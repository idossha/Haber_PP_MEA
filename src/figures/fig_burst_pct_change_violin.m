function fig_burst_pct_change_violin(study, varargin)
%FIG_BURST_PCT_CHANGE_VIOLIN Channel-wise burst-rate percent change violin.
%
%   FIG_BURST_PCT_CHANGE_VIOLIN(study) loads cached burst rates and spike
%   rates, applies the proof-of-concept burst-violin filters, computes
%   100*(t-b)/b per channel (only where baseline burst > 0), and writes
%   <study>_burst_pct_change_violin.png to output/fig{2,4}/panels/.
%
%   FIG_BURST_PCT_CHANGE_VIOLIN(study,                              ...
%        'pctMaxInclude',         1000,                              ...
%        'excludeSilencedPct',    false,                             ...
%        'ignoreSilentChannels',  true,                              ...
%        'minRateThreshold',      0,                                 ...
%        'pctYlim',               [],                                ...
%        'channels',              1:64)
%   overrides individual options.
%
% Channel filter (matches POC verbatim):
%   - Drop NaNs across (burst baseline, burst treatment, spike baseline,
%     spike treatment).
%   - Apply ignoreSilentChannels / minRateThreshold on burst rates.
%   - Drop channels with both spike baseline AND treatment <= 0 or both 0.
%   - Keep only channels with baseline burst > 0.
%   - excludeSilencedPct + pctMaxInclude as in the spike violin.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',              cfg.channels.default);
    addParameter(p, 'pctMaxInclude',         cfg.pct.max_include);
    addParameter(p, 'excludeSilencedPct',    cfg.pct.exclude_silenced);
    addParameter(p, 'ignoreSilentChannels',  true);
    addParameter(p, 'minRateThreshold',      cfg.silent.min_rate_burst);
    addParameter(p, 'pctYlim',               cfg.pct.pct_ylim);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);

    % We need both burst rates and spike rates aligned per pair, so load
    % the full cache struct for each half (load_pair_metric pools only one
    % metric at a time).
    pctAll = [];
    for k = 1:numel(pairs)
        [loadedB, loadedT] = load_pair_cache(pairs(k), cfg);

        burstB = loadedB.burstRates(:);
        burstT = loadedT.burstRates(:);
        spikeB = loadedB.spikeRates(:);
        spikeT = loadedT.spikeRates(:);
        chB    = loadedB.channelsUsed(:)';
        chT    = loadedT.channelsUsed(:)';

        [~, idxB] = ismember(opt.channels, chB);
        [~, idxT] = ismember(opt.channels, chT);
        validIdx  = (idxB > 0) & (idxT > 0);
        if ~any(validIdx); continue; end

        b  = burstB(idxB(validIdx));
        t  = burstT(idxT(validIdx));
        sB = spikeB(idxB(validIdx));
        sT = spikeT(idxT(validIdx));

        ok = ~isnan(b) & ~isnan(t) & ~isnan(sB) & ~isnan(sT);
        b  = b(ok); t = t(ok); sB = sB(ok); sT = sT(ok);

        if opt.ignoreSilentChannels
            if opt.minRateThreshold == 0
                nonSilent = ~(b == 0 & t == 0);
            else
                nonSilent = ~(b < opt.minRateThreshold & t < opt.minRateThreshold);
            end
            b = b(nonSilent); t = t(nonSilent); sB = sB(nonSilent); sT = sT(nonSilent);
        end
        if isempty(b); continue; end

        keepSpike = ~((sB < 0 & sT < 0) | (sB == 0 & sT == 0));
        b = b(keepSpike); t = t(keepSpike);

        use = b > 0;
        b = b(use); t = t(use);

        pctPair = 100 * (t - b) ./ b;
        pctAll = [pctAll; pctPair]; %#ok<AGROW>
    end

    if isempty(pctAll)
        error('fig_burst_pct_change_violin:NoData', 'No observations remain.');
    end

    if opt.excludeSilencedPct
        keep = abs(pctAll + 100) >= 1e-6;
        pctAll = pctAll(keep);
    end
    if ~isempty(opt.pctMaxInclude)
        keep = pctAll <= opt.pctMaxInclude;
        pctAll = pctAll(keep);
    end
    if isempty(pctAll)
        error('fig_burst_pct_change_violin:NoData', 'Filters left no observations.');
    end

    panelDir = output_path(cfg, study, 'rates', 'panels');
    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    colors  = paired_plot_colors(study);
    outFile = fullfile(panelDir, sprintf('%s_burst_pct_change_violin.png', study));

    plot_pct_change_violin(pctAll, colors, labels.treatment, ...
        'Percent change in burst rate (%)', opt.pctYlim, outFile);

    fprintf('fig_burst_pct_change_violin(%s): saved %s\n', study, outFile);

    stats = struct( ...
        'study',          study, ...
        'metric',         'burst_rate_pct_change', ...
        'n_observations', numel(pctAll), ...
        'median_pct',     median(pctAll), ...
        'iqr_q1',         prctile(pctAll, 25), ...
        'iqr_q3',         prctile(pctAll, 75), ...
        'mean_pct',       mean(pctAll), ...
        'pct_above_zero', 100 * mean(pctAll > 0), ...
        'pct_below_zero', 100 * mean(pctAll < 0), ...
        'figure_file',    outFile);
    export_figure_stats(stats, fullfile(statsDir, ...
        sprintf('%s_burst_pct_change_violin_stats', study)));
end
