# Haber_PP_MEA — analysis source

MATLAB pipeline for the DOI / ketanserin psychoplastogen MEA paper. Raw TDT
recordings are filtered, spike-detected, and burst-detected once into a per-
dataset cache. Every figure and analysis script reads from cache.

## Layout

```
src/
├── config/
│   └── project_config.m          # Single source of truth: paths, filters,
│                                 # spike/burst params, dataset lists, channels.
├── utils/
│   ├── get_project_root.m        # mfilename-based repo root resolver.
│   ├── build_notch_filter.m      # 60 Hz notch stack with harmonics up to fmax.
│   ├── build_bandpass_filter.m   # 4th-order Butterworth bandpass.
│   ├── load_and_filter.m         # TDTbin2mat + notch + bandpass for one channel.
│   ├── detect_spikes.m           # TDTthresh wrapper -> sorted spike times.
│   ├── detect_bursts.m           # ISI-threshold burst state machine.
│   ├── cache_filename.m          # dataset name -> cache_<safe>.mat
│   ├── dataset_path.m            # dataset name -> absolute folder path.
│   ├── mea60_layout.m            # 60-electrode MEA geometry (channel -> grid).
│   ├── pair_datasets.m           # Pair baselines and treatments by "#N" suffix.
│   ├── get_pairs_and_labels.m    # Convenience: pairs + axis labels per study.
│   ├── load_cache.m              # Load one dataset cache (schema-validated).
│   ├── load_pair_cache.m         # Load both halves of a pair (two structs).
│   ├── load_pair_metric.m        # Pool a per-channel summary metric across pairs.
│   ├── load_pair_spikes.m        # Pool spike-time cells for a pair (connectivity).
│   ├── nice_tick_step.m          # 1/2/5 * 10^n tick spacing.
│   ├── plot_mean_iqr_marker.m    # Tukey IQR bar + mean dot.
│   ├── plot_paired_lines.m       # Paired-line + half-violin figure.
│   ├── plot_paired_boxplot.m     # Two-group boxplot with Tukey-mean overlay.
│   ├── plot_pct_change_violin.m  # Vertical violin + boxplot for % change.
│   ├── paired_plot_colors.m      # Per-study color palette.
│   └── save_figure.m             # 300 dpi PNG with print() fallback.
├── pipeline/
│   └── preprocess_and_save.m     # Build per-dataset caches in cfg.paths.cache.
├── figures/
│   ├── fig_spike_rate.m              # Paired-line + boxplot of spikes/min.
│   ├── fig_burst_rate.m              # Paired-line + boxplot of bursts/min.
│   ├── fig_spike_pct_change_violin.m # Channel-wise spike % change violin.
│   ├── fig_burst_pct_change_violin.m # Channel-wise burst % change violin.
│   ├── fig_rate_change_bar.m         # Silenced/decreased/increased boxplots
│   │                                 # (metric = 'spike' | 'burst').
│   ├── fig_filtered_traces.m         # Side-by-side filtered voltage traces.
│   ├── fig_connectivity_exemplar.m   # Xcorrelograms + adjacency + network graph.
│   ├── fig_connectivity_summary.m    # Paired network-metric plots across pairs.
│   └── fig_stats_bootstrap.m         # Hierarchical bootstrap on % change.
├── analysis/
│   ├── connectivity_xcorr.m      # Pairwise binned spike-train cross-correlation.
│   ├── run_connectivity.m        # Driver: load pair spikes -> xcorr -> metrics.
│   ├── plot_network_on_mea.m     # Draw nodes + edges over the 60MEA geometry.
│   ├── topographical_map.m       # 8x8 MEA grid heatmap of any per-channel metric.
│   └── network_metrics.m         # Density, clustering, path length, global
│                                 # efficiency, modularity, small-worldness σ.
├── archive/
│   ├── oldpipeline_scripts/      # Original POC pipeline (read-only reference).
│   ├── other_scripts/            # Exploratory scripts (correlation, raster, etc.).
│   └── figure_scripts/           # Original POC figure scripts.
└── README.md
```

Nothing in `src/` (outside `archive/`) hardcodes a path or dataset name —
everything routes through `config/project_config.m`.

## Workflow

1. **Edit `src/config/project_config.m`** if you need to:
   - relocate the project (it auto-resolves via `get_project_root.m`),
   - change a filter / spike-detection / burst-detection parameter,
   - add a new dataset (append the folder name to the matching
     `cfg.datasets.<study>.baseline` and `treatment` cell arrays).

2. **Build the cache once** (run from any cwd):
   ```matlab
   addpath(genpath('src'));
   preprocess_and_save();                       % everything
   preprocess_and_save('study', 'doi');         % DOI only
   preprocess_and_save('overwrite', true);      % rebuild
   ```
   Caches land in `cfg.paths.cache` (`cache/cache_<dataset>.mat` by default).

3. **Render any figure** (DOI or ketanserin study):
   ```matlab
   addpath(genpath('src'));
   fig_spike_rate('doi');
   fig_spike_rate('ket');
   fig_burst_rate('doi');
   fig_spike_pct_change_violin('doi');
   fig_burst_pct_change_violin('doi');
   fig_rate_change_bar('spike', 'doi');
   fig_rate_change_bar('burst', 'doi');
   fig_filtered_traces('IdoControl-230914-130200_#1', ...
                       'IdoDOI-230914-142502_#1', [35 37 23 24 25]);
   stats = fig_stats_bootstrap('doi');
   ```
   PNGs are written to `cfg.paths.figures_out` (default `figures_out/`).

4. **Analysis modules**:
   ```matlab
   addpath(genpath('src'));
   cfg = project_config();

   % Connectivity for every pair in a study (uses cached spike times,
   % never re-runs preprocessing):
   results = run_connectivity('doi');                % all DOI pairs
   results = run_connectivity('doi', 'save', true);  % + persist per-pair
                                                     %   analysis caches

   % Or, for a one-off pair, load spikes directly and call the primitive:
   [pairs, ~]       = get_pairs_and_labels(cfg, 'doi');
   [bSp, tSp, bMt]  = load_pair_spikes(pairs(1), cfg.channels.default, cfg);
   conn             = connectivity_xcorr(bSp, 'durationSec', bMt.durationSec);
   metrics          = network_metrics(conn.adjacency);

   % Topographical map (per-channel metric -> 8x8 MEA grid):
   cache = load_cache(pairs(1).baseline, cfg);
   topographical_map(cache.spikeRates, 'cbarLabel', 'spikes/min', ...
       'outFile', fullfile('figures_out', 'topo_baseline.png'));
   ```

   See `docs/CACHE_SCHEMA.md` for the full on-disk format and a list of
   the derivation one-liners (ISI distributions, in-burst fraction, etc.).

## Adding a new dataset

1. Drop the TDT block folder under `data/DOI/` or `data/ketanserin/`. The
   folder name should end in `_#<integer>` so `pair_datasets()` can match
   it to its baseline / treatment counterpart.
2. Append the folder name to the matching cell arrays in
   `cfg.datasets.<study>.baseline` and `cfg.datasets.<study>.treatment` in
   `src/config/project_config.m` (preserve pair order).
3. Re-run `pipeline.preprocess_and_save()` (existing caches are skipped by
   default; the new dataset is processed once).
4. Re-run any figure / analysis script.

## Dependencies

- MATLAB R2020a+ (uses `exportgraphics`, `histcounts`, `xcorr`, `filtfilt`,
  `designNotchPeakIIR`, `boxplot`, `ksdensity`).
- The TDT MATLAB SDK lives at `<repo>/TDTMatlabSDK/`. The pipeline calls
  `addpath(genpath(cfg.paths.tdt_sdk))` automatically.
- Statistics and Machine Learning Toolbox (`prctile`, `chi2cdf`, `ksdensity`).
- Signal Processing Toolbox (`butter`, `filtfilt`, `designNotchPeakIIR`).

## Canonical processing parameters (see `project_config.m`)

| Parameter           | Value |
|---------------------|-------|
| Notch fundamental   | 60 Hz |
| Notch harmonics     | up to 780 Hz (zero-phase, Q = 35) |
| Bandpass            | 300-2500 Hz, 4th-order Butterworth, zero-phase |
| Spike detection     | `TDTthresh` MODE=auto, POLARITY=-1, STD=6.5, TAU=5 |
| Burst detection     | ISI <= 100 ms starts, ISI > 200 ms ends, min 3 spikes |
| Channels            | 1:60 (59 recording + 1 iR reference) |

## Archive

`src/archive/` preserves the original proof-of-concept code (do not modify).
The connectivity inspiration for `analysis/connectivity_xcorr.m` came from
`archive/other_scripts/Ido_Ilhan__Correlation_automatic*.m` and
`archive/other_scripts/OnlyCross_data_adaptive_zero_delay.m`. The
topographical channel-to-grid mapping came from
`archive/other_scripts/CartesianPlot_histo_Rasters_Dataset.m`. The burst
detector matches `archive/oldpipeline_scripts/burst_stats.m` verbatim.
