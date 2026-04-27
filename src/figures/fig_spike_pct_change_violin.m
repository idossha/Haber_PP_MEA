function fig_spike_pct_change_violin(study, varargin)
%FIG_SPIKE_PCT_CHANGE_VIOLIN Channel-wise spike-rate percent change violin.
%
%   FIG_SPIKE_PCT_CHANGE_VIOLIN(study) loads cached spike rates for the
%   requested study, computes 100*(t-b)/b per channel using the same
%   filters as the proof-of-concept figures_spike_pct_change_violin.m, and
%   writes <study>_spike_pct_change_violin.png to output/fig{2,4}/panels/.
%
%   FIG_SPIKE_PCT_CHANGE_VIOLIN(study,                          ...
%        'pctMaxInclude',     1000,                              ...
%        'excludeSilencedPct', false,                            ...
%        'pctYlim',           [],                                ...
%        'channels',          1:64)
%   overrides individual options.
%
% Channel filter (matches POC verbatim):
%   - Drop channels with both baseline and treatment spike rate == 0.
%   - Inf percent change (b == 0, t > 0) is capped at pctMaxInclude (or 1e6
%     when pctMaxInclude is empty).
%   - When excludeSilencedPct is true, drop -100% observations.
%   - When pctMaxInclude is non-empty, drop pct > pctMaxInclude.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',           cfg.channels.default);
    addParameter(p, 'pctMaxInclude',      cfg.pct.max_include);
    addParameter(p, 'excludeSilencedPct', cfg.pct.exclude_silenced);
    addParameter(p, 'pctYlim',            cfg.pct.pct_ylim);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    [bRates, tRates] = load_pair_metric(pairs, opt.channels, 'spikeRates', cfg);

    if isempty(bRates)
        error('fig_spike_pct_change_violin:NoData', 'No aligned channels.');
    end

    % Drop both-zero pairs (keep b==0,t>0 as Inf percent change).
    keep   = ~(bRates == 0 & tRates == 0);
    bRates = bRates(keep);
    tRates = tRates(keep);
    pctAll = 100 * (tRates - bRates) ./ bRates;

    infMask = isinf(pctAll);
    if any(infMask)
        if ~isempty(opt.pctMaxInclude)
            pctAll(infMask) = opt.pctMaxInclude;
        else
            pctAll(infMask) = 1e6;
        end
    end

    if opt.excludeSilencedPct
        nonSilent = abs(pctAll + 100) >= 1e-6;
        nDrop = numel(pctAll) - sum(nonSilent);
        pctAll = pctAll(nonSilent);
        fprintf('Dropped %d silenced (-100%%) observations.\n', nDrop);
    end

    if ~isempty(opt.pctMaxInclude)
        keepHi = pctAll <= opt.pctMaxInclude;
        nDropHi = sum(~keepHi);
        pctAll = pctAll(keepHi);
        fprintf('Dropped %d observations above pctMaxInclude=%g.\n', nDropHi, opt.pctMaxInclude);
    end

    if isempty(pctAll)
        error('fig_spike_pct_change_violin:NoData', 'No observations remain after filters.');
    end

    panelDir = output_path(cfg, study, 'rates', '');
    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    colors  = paired_plot_colors(study);
    outFile = fullfile(panelDir, sprintf('%s_spike_pct_change_violin.png', study));

    plot_pct_change_violin(pctAll, colors, labels.treatment, ...
        'Percent change in firing rate (%)', opt.pctYlim, outFile);

    fprintf('fig_spike_pct_change_violin(%s): saved %s\n', study, outFile);

    stats = struct( ...
        'study',          study, ...
        'metric',         'spike_rate_pct_change', ...
        'n_observations', numel(pctAll), ...
        'median_pct',     median(pctAll), ...
        'iqr_q1',         prctile(pctAll, 25), ...
        'iqr_q3',         prctile(pctAll, 75), ...
        'mean_pct',       mean(pctAll), ...
        'pct_above_zero', 100 * mean(pctAll > 0), ...
        'pct_below_zero', 100 * mean(pctAll < 0), ...
        'figure_file',    outFile);
    export_figure_stats(stats, fullfile(statsDir, ...
        sprintf('%s_spike_pct_change_violin_stats', study)));
end
