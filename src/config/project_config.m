function cfg = project_config()
%PROJECT_CONFIG Single source of truth for paths, filters, and dataset lists.
%
%   cfg = PROJECT_CONFIG() returns a struct with every parameter, dataset
%   list, and default used by the preprocessing pipeline, figure scripts,
%   and analysis functions. Edit this file to change defaults; nothing
%   else in src/ should hardcode any of these values.
%
%   EVERY non-obvious default is annotated with the primary literature
%   citation that supports it. For the full write-up and sensitivity
%   analysis, see:
%       docs/DECISION_VARIABLES.md
%       docs/METRICS.md
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   cfg  -  Struct with nested fields:
%             cfg.paths       - root, data, cache, figures_out, tdt_sdk
%             cfg.figures     - output format, DPI
%             cfg.filter      - notch + bandpass parameters
%             cfg.spike       - TDTthresh parameters
%             cfg.burst       - ISI-threshold burst parameters
%             cfg.silent      - silent-channel handling
%             cfg.outlier     - robust outlier filter
%             cfg.pct         - percent-change behaviour
%             cfg.connectivity- xcorr + network_metrics defaults
%             cfg.stats       - bootstrap, FDR, null ensemble
%             cfg.tdt         - TDT-specific constants
%             cfg.datasets    - study -> pair lists
%             cfg.labels      - axis labels per study
%             cfg.channels    - default channel list
%
% Notes:
%   - Dataset names are folder names (without trailing slash) inside the
%     corresponding cfg.paths.data_* directory.
%   - Pairing convention: cfg.datasets.<study>.baseline{k} is the matched
%     control for cfg.datasets.<study>.treatment{k}; pair_datasets() will
%     verify pairing by the trailing "#N" suffix.

% =========================================================================
% Paths
% =========================================================================
root                  = get_project_root();
cfg.paths.root        = root;
cfg.paths.data_doi    = fullfile(root, 'data', 'DOI');
cfg.paths.data_ket    = fullfile(root, 'data', 'ketanserin');
cfg.paths.cache       = fullfile(root, 'cache');
cfg.paths.figures_out = fullfile(root, 'figures', 'panels');
cfg.paths.tdt_sdk     = fullfile(root, 'TDTMatlabSDK');

% =========================================================================
% Figure output
% =========================================================================
% 'pdf+png'   writes both (canonical; .pdf for Inkscape editing, .png for preview)
% 'pdf'       PDF only
% 'png'       PNG only
cfg.figures.format    = {'pdf','png'};
cfg.figures.dpi       = 600;    % raster resolution for PNG + rasterised PDF elements

% =========================================================================
% TDT acquisition
% =========================================================================
cfg.tdt.store = 'Wav1';

% =========================================================================
% Signal-processing filters
% =========================================================================
% Sampling rate is read from each dataset; this is the TDT default and is
% never used as authoritative.
cfg.filter.fs_default = 24414.0625;

% 60 Hz notch stack (zero-phase, harmonics up to notch_max).
%   f0 = 60                     [standard for North American mains]
%   notch_max = 780 Hz (13 h.)  [Chiappalone 2006 / Martinoia 2005 — high-Q
%                                notches needed because the spike band
%                                starts at 300 Hz]
%   Q = 35                      [high Q to avoid widening the stop-band
%                                into the spike band]
cfg.filter.notch_f0  = 60;
cfg.filter.notch_max = 780;
cfg.filter.notch_Q   = 35;

% Bandpass for spike-band signal (zero-phase 4th-order Butterworth).
%   300-2500 Hz is the Chiappalone 2006 / Martinoia 2005 convention.
%   4th order is standard for MEA spike-band filtering.
cfg.filter.bp_low   = 300;
cfg.filter.bp_high  = 2500;
cfg.filter.bp_order = 4;

% =========================================================================
% Spike detection (TDTthresh)
% =========================================================================
% MODE='auto'          robust noise estimation per channel
% POLARITY=-1          negative-going spikes
% STD=5.0              threshold = 5.0 * robust SD of noise
%                      PRIMARY: Schroeter 2015 p.5460, Kapucu 2016 p.3,
%                      Hernandes 2024 p.2 all use 5 SD. 5 is the modal
%                      literature default. 6.5 is retained as a
%                      sensitivity supplement.
% TAU=5                seconds of MOVING-WINDOW used by TDTthresh to
%                      estimate the per-channel baseline noise floor
%                      (NOT a refractory period; TDTthresh documentation
%                      line 35-36: "defines size of moving window for
%                      'auto' thresholding mode, in seconds"). 5 s is
%                      the TDT default and matches the POC.
% refractory_ms=1.0    Post-detection refractory window enforced by
%                      detect_spikes.m (Brofiga 2023 p.4: 1 ms;
%                      Garofalo 2009 p.3: 1 ms). Spikes within this
%                      window of the previous detected spike are
%                      removed. Set to 0 to disable.
cfg.spike.mode          = 'auto';
cfg.spike.polarity      = -1;
cfg.spike.std_factor    = 5.0;
cfg.spike.tau           = 5;     % seconds (TDTthresh noise window)
cfg.spike.refractory_ms = 1.0;   % post-detection refractory (ms)

% =========================================================================
% Burst detection (ISI threshold; matches oldpipeline/burst_stats.m)
% =========================================================================
% Chiappalone 2006 uses ISI_max = 100 ms; Pasquale 2010 builds on that.
% min_spikes: Brofiga 2023 p.4 uses 5 explicitly; Chiappalone 2006 uses
% 5-10 depending on variant; Bakkum 2014 uses 10 for ISI_N. 3 is on the
% inclusive edge of the literature. Primary = 5; keep 3 as sensitivity
% supplement for reviewer response.
cfg.burst.isi_max_ms = 100;   % ISI <= this starts or continues a burst
cfg.burst.isi_end_ms = 200;   % ISI >  this ends a burst
cfg.burst.min_spikes = 5;     % primary (Brofiga 2023)

% =========================================================================
% Silent-channel handling
% =========================================================================
% Applied to spike rate (spikes/min) before paired analyses.
%   mode = 'baseline_min'   PRIMARY. Drop channels where baseline spike
%                           rate < min_rate_spike (ignores treatment).
%                           Mossink 2021 p.2187: "MFR > 0.1 spike/s".
%                           Brofiga 2023 p.5: "if MFR < 0.1 spikes/s,
%                           channel was discarded". This is the single
%                           most-citable silent-channel default in the
%                           MEA pharmacology literature.
%   mode = 'both_zero'      sensitivity supplement: drop only when BOTH
%                           baseline and treatment rate == 0. Permissive
%                           rule used by our POC / Varley 2024.
%   mode = 'either_zero'    most aggressive.
%   min_rate_spike = 6.0    = 0.1 spikes/s * 60 s/min (Mossink 2021).
%   min_rate_burst = 0.0    drop only (0, 0) burst-pair rows.
cfg.silent.mode            = 'baseline_min';
cfg.silent.min_rate_spike  = 6.0;     % Mossink 2021 / Brofiga 2023
cfg.silent.min_rate_burst  = 0.0;

% =========================================================================
% Outlier filter (robust_outlier_filter.m)
% =========================================================================
% IMPORTANT: after surveying the MEA pharmacology literature the research
% agent concluded: **no mainstream MEA paper removes high-end rate
% outliers**. Mossink 2021, Brofiga 2023, Varley 2024, and Olson 2024 all
% apply ONLY a low-end silent-channel cut and then run non-parametric
% stats; Brofiga 2023 p.5 and Varley 2024 log-transform asymmetric
% distributions for display instead of dropping points. The Tukey IQR
% rule was designed for roughly symmetric distributions and is
% inappropriate for heavy-tailed spike-rate data.
%
% DEFAULT: do NOT drop outliers for the primary analysis. Violin figures
% use Y-axis clipping (cfg.pct.violin_ylim_primary) instead. If a figure
% IS filtered (e.g. for a supplementary "cleaner" view) use the
% "extreme" Tukey 3.0 x IQR fence, not the 1.5 x "mild" fence.
%
%   mode = 'none'              PRIMARY (literature convention)
%   mode = 'tukey'              supplementary: iqr_factor = 3.0
%   mode = 'percentile'         supplementary: drop > upper_percentile
%   mode = 'mean_multiplier'    legacy POC mode
cfg.outlier.mode             = 'none';
cfg.outlier.iqr_factor       = 3.0;   % "extreme" if ever used
cfg.outlier.upper_percentile = 99;
cfg.outlier.mean_multiplier  = 15;

% =========================================================================
% Percent change and violins
% =========================================================================
% max_include          cap (%) at which b=0, t>0 observations enter the
%                      pooled distribution. 1000 % is the standard cap.
% exclude_silenced     drop observations with pct ~ -100% (complete
%                      silencing) from violin plots. Default false:
%                      silencing is valid biological data.
% violin_ylim_primary  visual y-range for MAIN-TEXT figures. Replaces
%                      the need for outlier removal; the underlying
%                      pooled distribution is not truncated, just
%                      clipped for display. -100 to +500 is the
%                      convention (Brofiga 2023 log-transform analogue,
%                      Mossink 2021 rate-change plots).
% violin_ylim_supplement  wider y-range for supplementary cap-sensitivity
%                      figures.
cfg.pct.max_include            = 1000;
cfg.pct.exclude_silenced       = false;
cfg.pct.pct_ylim               = [-100 300];    % primary (tightened for publication)
cfg.pct.violin_ylim_supplement = [-100 1000];   % supplementary view

% =========================================================================
% Connectivity (connectivity_xcorr.m + network_metrics.m)
% =========================================================================
% bin_ms           1        Chiappalone 2006, Garofalo 2009, Pastore 2016
% max_lag_ms       100      primary (Chiappalone 2006, Schroeter 2015).
%                           Also report a ±50 ms "monosynaptic" supplement
%                           (Schroeter 2015 p.5461 uses ±15 ms).
% normalization    'pearson' z-scored binned spike trains (primary).
%                           STTC supplement added separately (Cutts &
%                           Eglen 2014) because Pearson-on-counts fails
%                           the rate-independence test (Cutts 2014,
%                           Table 3).
% edge_density     0.10     PRIMARY value for main-text figures (Downes
%                           2012 uses 10 %). Supplement sweeps
%                           [0.05 0.10 0.15 0.20 0.25] per Rubinov &
%                           Sporns 2010 recommendation.
% edge_threshold   []       absolute threshold overrides data-driven
%                           percentile from edge_density; [] = use density
% density_sweep    ...       reported in supplementary figure
% null_N           1000     null ensemble size for small-worldness sigma
%                           (Humphries & Gurney 2008 p.2 uses 1000
%                           Erdos-Renyi graphs for 99 % CIs). For faster
%                           iteration during development use 100.
% null_seed        1        RNG seed for reproducibility
cfg.connectivity.bin_ms         = 1;
cfg.connectivity.max_lag_ms     = 100;
cfg.connectivity.max_lag_ms_supp = 50;   % monosynaptic supplement
cfg.connectivity.normalization  = 'pearson';
cfg.connectivity.edge_density   = 0.10;   % Downes 2012 primary
cfg.connectivity.edge_threshold = [];
cfg.connectivity.density_sweep  = [0.05 0.10 0.15 0.20 0.25];
cfg.connectivity.null_N         = 100;    % 100 primary, 1000 for final submission
cfg.connectivity.null_seed      = 1;

% =========================================================================
% Statistics
% =========================================================================
% bootstrap iterations: 10000 is the Saravanan 2020 recommendation for
% median stability; 1000 is the minimum defensible; 10000 is standard.
cfg.stats.n_bootstrap          = 10000;
cfg.stats.bootstrap_seed       = 1;
cfg.stats.fdr_q                = 0.05;   % Benjamini-Hochberg FDR target
cfg.stats.min_n_wilcoxon_exact = 20;     % below this, use exact signrank

% =========================================================================
% Avalanches (Pasquale 2008 + Massobrio 2015)
% =========================================================================
% Neither Pasquale 2008 nor Massobrio 2015 defines a branching ratio sigma
% formula or a Deviation-from-Criticality Coefficient (DCC) — see
% Tracks/Active/phase1_decisions.md §1. This block anchors ONLY to values
% that are explicitly backed by a page citation in those two PDFs. sigma
% and DCC are deliberately NOT reported; the Sethna scaling residual
% (Massobrio 2015 p.9 Eq. 1-2) is the closest-supported criticality-
% consistency scalar and is computed instead.
cfg.avalanches.bin_selection            = 'mean_iei';                     % Massobrio 2015 p.2, p.8
cfg.avalanches.bin_sweep_ms             = [0.2 0.4 0.6 0.8 1 2 4 8 16];   % Pasquale 2008 p.1357
cfg.avalanches.active_threshold_spikes  = 1;                               % Pasquale 2008 p.1357; Massobrio 2015 p.2
cfg.avalanches.size_definition          = 'unique_electrodes';             % Pasquale 2008 p.1357 (defn 2), p.1361
cfg.avalanches.drop_boundary_avalanches = true;                            % Pasquale 2008 p.1357 flanking-silent-bins rule
cfg.avalanches.alpha_size_target        = -1.5;                            % Pasquale 2008 p.1358, p.1361; Massobrio 2015 p.2
cfg.avalanches.beta_lifetime_target     = -2.0;                            % Pasquale 2008 p.1358, p.1361
cfg.avalanches.alpha_tolerance          = 0.15;                            % Pasquale 2008 p.1361 empirical SD 0.09-0.13
cfg.avalanches.tail_mass_frac           = 0.80;                            % Pasquale 2008 Fig. 2/6 tail-upturn rule
cfg.avalanches.ls_drop_unit_bin         = true;                            % Pasquale 2008 p.1357
cfg.avalanches.ls_drop_below_pct_of_max = 0.01;                            % Pasquale 2008 p.1357
cfg.avalanches.mle_enabled              = true;                            % Massobrio 2015 p.14
cfg.avalanches.mle_alternatives         = {'exponential','truncated_power_law','lognormal'}; % Massobrio 2015 p.14
cfg.avalanches.ks_pvalue_threshold      = 0.10;                            % Massobrio 2015 p.5
cfg.avalanches.compute_sethna_residual  = true;                            % Massobrio 2015 p.9 Eq. 1-2

% =========================================================================
% Transfer entropy (Ito 2011)
% =========================================================================
cfg.te.bin_ms         = 1;      % Ito 2011 p.4 (canonical D1TE 1 ms bin)
cfg.te.delay_bins     = 1;      % Ito 2011 p.4 (D1TE = delay-1)
cfg.te.min_rate_hz    = 5;      % Ito 2011 active-channel eligibility
cfg.te.top_edges_frac = 0.10;   % top-decile summary; tmp/NEW_ANGLES_REPORT.md §1

% =========================================================================
% Shannon entropy (Varley 2024)
% =========================================================================
cfg.entropy.bin_ms      = 50;   % Varley 2024 §2.2 convention
cfg.entropy.clip_counts = 3;    % Varley 2024 §2.2 convention (alphabet {0..3}, H_max = 2 bits)

% =========================================================================
% Burst-onset coincidence (Brofiga 2023 + Chiappalone 2006)
% =========================================================================
cfg.burst_sync.window_ms       = 50;                % Chiappalone 2006 (burst propagation <100 ms)
cfg.burst_sync.window_sweep_ms = [25 50 100 200];   % SI §S4 sweep

% =========================================================================
% Channels
% =========================================================================
cfg.channels.default = 1:64;

% =========================================================================
% Dataset lists
% =========================================================================
% IMPORTANT: include ALL datasets (no exclusions). Pairs are matched by
% trailing "#N" suffix. To add a new dataset, append it to the matching
% baseline/treatment cell array below, in pair order.

% DOI study (data/DOI/) - 6 paired recordings
cfg.datasets.doi.baseline = { ...
    'IdoControl-230914-130200_#1', ...
    'IdoControl-230914-131548_#2', ...
    'IdoControl-230914-132855_#3', ...
    'IdoControl-230914-154022_#4', ...
    'IdoControl-230914-155318_#5', ...
    'IdoControl-230914-160601_#6'};
cfg.datasets.doi.treatment = { ...
    'IdoDOI-230914-142502_#1', ...
    'IdoDOI-230914-143740_#2', ...
    'IdoDOI-230914-144945_#3', ...
    'IdoDOI-230914-161838_#4', ...
    'IdoDOI-230914-163140_#5', ...
    'IdoDOI-230914-164512_#6'};

% Ketanserin study (data/ketanserin/) - 3 paired recordings
cfg.datasets.ket.baseline = { ...
    'IdoControl-Kt-231101-144046_#1', ...
    'IdoControl-Kt-231101-145330_#2', ...
    'IdoControl-Kt-231101-150512_#3'};
cfg.datasets.ket.treatment = { ...
    'IdoKetanserin-231101-133725_#1', ...
    'IdoKetanserin-231101-134946_#2', ...
    'IdoKetanserin-231101-140237_#3'};

% Convenience: per-study labels for axes/legends.
cfg.labels.doi.baseline  = 'Baseline';
cfg.labels.doi.treatment = 'DOI';
cfg.labels.ket.baseline  = 'Baseline';
cfg.labels.ket.treatment = 'Ketanserin';

end
