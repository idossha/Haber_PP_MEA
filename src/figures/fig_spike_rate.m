function fig_spike_rate(study, varargin)
%FIG_SPIKE_RATE Paired-line + boxplot of spike rate (spikes/min) per channel.
%
%   FIG_SPIKE_RATE(study) loads the cached spike rates for the requested
%   study ('doi' or 'ket'), pools all baseline/treatment pairs, applies the
%   canonical channel filters (silent-channel exclusion, mean-multiplier
%   outlier filter), and writes two figures to output/fig{2,4}/panels/:
%       <study>_spike_rate_paired.png    -  paired-line + half-violin
%       <study>_spike_rate_boxplot.png   -  paired boxplot
%
%   FIG_SPIKE_RATE(study, 'minRateThreshold', 1.0,                ...
%                         'ignoreSilentChannels', true,            ...
%                         'excludeMeanMultiplierOutliers', true,   ...
%                         'meanRateOutlierMultiplier', 15,         ...
%                         'channels', 1:64)
%   overrides individual options. Defaults match the proof-of-concept
%   src/figure_scripts/figures_spike_rate.m exactly.
%
% INPUTS:
%   study  -  'doi' or 'ket'.
%
% OUTPUTS:
%   PNG files written to output/fig{2,4}/panels/.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',                       cfg.channels.default);
    addParameter(p, 'ignoreSilentChannels',           true);
    addParameter(p, 'silentMode',                     cfg.silent.mode, ...
        @(s) any(strcmpi(s, {'baseline_min','both_zero','either_zero'})));
    addParameter(p, 'minRateThreshold',               cfg.silent.min_rate_spike);
    addParameter(p, 'excludeMeanMultiplierOutliers',  false);
    addParameter(p, 'meanRateOutlierMultiplier',      cfg.outlier.mean_multiplier);
    addParameter(p, 'outlierMode',                    cfg.outlier.mode, ...
        @(s) any(strcmpi(s, {'none','tukey','percentile','mean_multiplier'})));
    addParameter(p, 'outlierUpperPercentile',         cfg.outlier.upper_percentile);
    addParameter(p, 'outlierIqrFactor',               cfg.outlier.iqr_factor);
    addParameter(p, 'yScale',                         'linear', ...
        @(s) any(strcmpi(s, {'linear','log'})));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    [ratesBaseline, ratesTreatment] = load_pair_metric(pairs, opt.channels, 'spikeRates', cfg);

    if isempty(ratesBaseline)
        error('fig_spike_rate:NoData', ...
            'No valid channels found across pairs after NaN drop.');
    end

    % --- Silent-channel filter -------------------------------------------
    % Primary literature default is 'baseline_min' at 0.1 spike/s (=
    % 6 spike/min): Mossink 2021 p.2187, Brofiga 2023 p.5. Kept the
    % older 'both_zero' / 'either_zero' rules as options for sensitivity
    % supplements and for backwards compat with existing caller scripts.
    nBefore = numel(ratesBaseline);
    if opt.ignoreSilentChannels
        switch lower(opt.silentMode)
            case 'baseline_min'
                nonSilent = ratesBaseline >= opt.minRateThreshold;
            case 'both_zero'
                nonSilent = ~(ratesBaseline < opt.minRateThreshold ...
                            & ratesTreatment < opt.minRateThreshold);
            case 'either_zero'
                nonSilent = ratesBaseline  >= opt.minRateThreshold ...
                          & ratesTreatment >= opt.minRateThreshold;
            otherwise
                nonSilent = true(size(ratesBaseline));
        end
        ratesBaseline  = ratesBaseline(nonSilent);
        ratesTreatment = ratesTreatment(nonSilent);
    end
    nSilentDropped = nBefore - numel(ratesBaseline);
    fprintf('Silent filter (%s, min=%.2f/min): dropped %d of %d.\n', ...
        opt.silentMode, opt.minRateThreshold, nSilentDropped, nBefore);
    % --- Outlier filtering ------------------------------------------------
    % The POC used a "mean multiplier" filter (drop anything > 15x mean).
    % The paper uses the more standard Tukey 1.5*IQR upper-fence filter
    % on pooled baseline+treatment rates (Chiappalone 2006, Mossink 2021).
    % Backward compat: the legacy mean_multiplier path still works if
    % either excludeMeanMultiplierOutliers is true or outlierMode is set
    % explicitly to 'mean_multiplier'.
    effectiveMode = opt.outlierMode;
    if opt.excludeMeanMultiplierOutliers && strcmpi(opt.outlierMode, 'tukey')
        % Legacy behaviour: legacy flag takes precedence over default mode
        % but not over an explicit user choice.
        effectiveMode = 'mean_multiplier';
    end
    [keepMask, outlierInfo] = robust_outlier_filter(ratesBaseline, ratesTreatment, ...
        'mode',            effectiveMode, ...
        'upperPercentile', opt.outlierUpperPercentile, ...
        'multiplier',      opt.meanRateOutlierMultiplier, ...
        'iqrFactor',       opt.outlierIqrFactor);
    ratesBaseline  = ratesBaseline(keepMask);
    ratesTreatment = ratesTreatment(keepMask);
    fprintf('Outlier filter (%s): dropped %d of %d channel observation(s) [cut=%.2f].\n', ...
        outlierInfo.mode, outlierInfo.nDropped, outlierInfo.nTotal, ...
        outlierInfo.upperCutBaseline);

    if isempty(ratesBaseline)
        error('fig_spike_rate:NoData', ...
            'All channels filtered out before plotting.');
    end

    nCh         = numel(ratesBaseline);
    pctIncrease = 100 * sum(ratesTreatment > ratesBaseline) / nCh;
    pctDecrease = 100 * sum(ratesTreatment < ratesBaseline) / nCh;

    panelDir = output_path(cfg, study, 'rates', 'panels');
    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    colors      = paired_plot_colors(study);
    pairedFile  = fullfile(panelDir, sprintf('%s_spike_rate_paired.png',  study));
    boxplotFile = fullfile(panelDir, sprintf('%s_spike_rate_boxplot.png', study));

    yLabelText = 'Firing rate (spikes min^{-1})';
    titleText  = sprintf('Firing rate: %s vs. %s', ...
                         labels.baseline, labels.treatment);

    plot_paired_lines(ratesBaseline, ratesTreatment, colors, ...
        '', ...
        pairedFile, ...
        'title',       titleText, ...
        'yLabel',      yLabelText, ...
        'xTickLabels', {labels.baseline, labels.treatment}, ...
        'yLimitMode',  'robust', ...
        'yScale',      opt.yScale, ...
        'yFloor',      5);

    plot_paired_boxplot(ratesBaseline, ratesTreatment, colors, labels, boxplotFile, ...
        'title',      titleText, ...
        'yLabel',     yLabelText, ...
        'yLimitMode', 'robust', ...
        'yScale',     opt.yScale);

    fprintf('fig_spike_rate(%s): n=%d, +%.1f%%, -%.1f%%\n', ...
        study, nCh, pctIncrease, pctDecrease);
    fprintf('  saved: %s\n  saved: %s\n', pairedFile, boxplotFile);

    % --- Numeric sidecar for paper writing ------------------------------
    psStats = paired_stats(ratesBaseline, ratesTreatment);
    stats = struct( ...
        'study',            study, ...
        'metric',           'spike_rate', ...
        'unit',             'spikes/min', ...
        'n_channels',       nCh, ...
        'pct_increased',    pctIncrease, ...
        'pct_decreased',    pctDecrease, ...
        'median_baseline',  psStats.medianBaseline, ...
        'median_treatment', psStats.medianTreatment, ...
        'median_delta',     psStats.medianDelta, ...
        'median_pct_change',psStats.medianPctChange, ...
        'ci_pct_change',    psStats.bootstrap.ciPctChange, ...
        'p_bootstrap',      psStats.bootstrap.pPctChange, ...
        'p_wilcoxon',       psStats.wilcoxon.p, ...
        'hedges_g_av',      psStats.hedgesGav, ...
        'figures_paired',   pairedFile, ...
        'figures_boxplot',  boxplotFile);
    export_figure_stats(stats, fullfile(statsDir, ...
        sprintf('%s_spike_rate_stats', study)));
end
