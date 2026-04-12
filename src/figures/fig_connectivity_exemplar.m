function fig_connectivity_exemplar(study, varargin)
%FIG_CONNECTIVITY_EXEMPLAR Multi-panel connectivity figure for one pair.
%
%   FIG_CONNECTIVITY_EXEMPLAR(study) picks one representative pair from
%   STUDY ('doi' or 'ket'), builds the panels of the connectivity figure,
%   and writes <study>_connectivity_exemplar.png to cfg.paths.figures_out.
%
%   Panel A  - two channel-pair cross-correlograms, baseline vs treatment.
%   Panel B  - adjacency heatmaps: baseline | treatment | delta.
%   Panel C  - network graph on the 60MEA geometry: baseline, treatment,
%              and the delta (threshold-based line colour).
%
%   FIG_CONNECTIVITY_EXEMPLAR(study, ...
%         'pairIndex',         3,       ...
%         'channels',          1:60,    ...
%         'binMs',             1,       ...
%         'maxLagMs',          100,     ...
%         'edgeThresholdPct',  20,      ...
%         'normalization',     'pearson')
%   overrides defaults. 'pairIndex' selects which baseline/treatment
%   pair to plot (default: middle pair).
%
% INPUTS:
%   study  -  'doi' | 'ket'
%
% OUTPUTS:
%   One PNG at <cfg.paths.figures_out>/<study>_connectivity_exemplar.png
%
% See also: CONNECTIVITY_XCORR, PLOT_NETWORK_ON_MEA, RUN_CONNECTIVITY.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pairIndex',        [],  @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'channels',         cfg.channels.default);
    addParameter(p, 'binMs',            cfg.connectivity.bin_ms,        @(x) isscalar(x) && x > 0);
    addParameter(p, 'maxLagMs',         cfg.connectivity.max_lag_ms,    @(x) isscalar(x) && x > 0);
    addParameter(p, 'normalization',    cfg.connectivity.normalization);
    addParameter(p, 'edgeThresholdPct', 100 * cfg.connectivity.edge_density, ...
        @(x) isscalar(x) && x >= 0 && x <= 100);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_connectivity_exemplar:NoPairs', 'No pairs for study %s.', study);
    end
    if isempty(opt.pairIndex)
        pairIdx = max(1, round(numel(pairs) / 2));
    else
        pairIdx = min(opt.pairIndex, numel(pairs));
    end
    pair = pairs(pairIdx);

    fprintf('fig_connectivity_exemplar(%s): pair %d of %d\n  baseline: %s\n  treatment: %s\n', ...
        study, pairIdx, numel(pairs), pair.baseline, pair.treatment);

    % --- Load cached spike times and per-channel rates ------------------
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, opt.channels, cfg);
    bCache = load_cache(pair.baseline,  cfg);
    tCache = load_cache(pair.treatment, cfg);

    % --- Run cross-correlation for both halves --------------------------
    xcorrArgs = { ...
        'binMs',         opt.binMs, ...
        'maxLagMs',      opt.maxLagMs, ...
        'normalization', opt.normalization};
    bRes = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, xcorrArgs{:});
    tRes = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, xcorrArgs{:});

    delta = tRes.adjacency - bRes.adjacency;

    % --- Per-channel spike rates aligned to opt.channels ----------------
    bRates = align_metric(bCache.spikeRates, bCache.channelsUsed, opt.channels);
    tRates = align_metric(tCache.spikeRates, tCache.channelsUsed, opt.channels);

    % --- Figure layout --------------------------------------------------
    colors = paired_plot_colors(study);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 1100]);

    tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Functional connectivity: %s vs %s (%s)', ...
          labels.baseline, labels.treatment, pair.baseline), ...
          'Interpreter', 'none');

    % ===== Panel A — two exemplar cross-correlograms ===================
    % Pick two pairs: (1) strongest baseline edge, (2) largest positive delta.
    [iTop, jTop]   = top_edge(bRes.adjacency);
    [iGain, jGain] = top_edge(delta);

    axA1 = nexttile(tl, 1);
    plot_crosscorrelogram_pair(axA1, bSpikes, tSpikes, ...
        iTop, jTop, bMeta.durationSec, tMeta.durationSec, ...
        opt, colors, labels, 'A1');
    title(axA1, sprintf('Strongest baseline edge  (ch %d \\leftrightarrow ch %d)', iTop, jTop));

    axA2 = nexttile(tl, 2);
    plot_crosscorrelogram_pair(axA2, bSpikes, tSpikes, ...
        iGain, jGain, bMeta.durationSec, tMeta.durationSec, ...
        opt, colors, labels, 'A2');
    title(axA2, sprintf('Largest \\Delta edge  (ch %d \\leftrightarrow ch %d)', iGain, jGain));

    axLegend = nexttile(tl, 3);
    axis(axLegend, 'off');
    text(axLegend, 0.05, 0.92, sprintf('study: %s', study), 'FontWeight', 'bold');
    text(axLegend, 0.05, 0.78, sprintf('pair #%d / %d', pairIdx, numel(pairs)));
    text(axLegend, 0.05, 0.66, sprintf('bin = %g ms', opt.binMs));
    text(axLegend, 0.05, 0.54, sprintf('max lag = %g ms', opt.maxLagMs));
    text(axLegend, 0.05, 0.42, sprintf('norm: %s', opt.normalization));
    text(axLegend, 0.05, 0.30, sprintf('edge threshold: top %.0f%%', opt.edgeThresholdPct));

    % ===== Panel B — adjacency heatmaps =================================
    axB1 = nexttile(tl, 4);
    plot_adjacency_heatmap(axB1, bRes.adjacency, 'Baseline', false);

    axB2 = nexttile(tl, 5);
    plot_adjacency_heatmap(axB2, tRes.adjacency, labels.treatment, false);

    axB3 = nexttile(tl, 6);
    plot_adjacency_heatmap(axB3, delta, '\Delta (treatment - baseline)', true);

    % ===== Panel C — network on MEA geometry ============================
    axC1 = nexttile(tl, 7);
    plot_network_on_mea(bRes.adjacency, ...
        'parent',           axC1, ...
        'nodeMetric',       bRates, ...
        'edgeMode',         'weight', ...
        'edgeThresholdPct', opt.edgeThresholdPct, ...
        'title',            sprintf('%s network (nodes = spikes/min)', labels.baseline));

    axC2 = nexttile(tl, 8);
    plot_network_on_mea(tRes.adjacency, ...
        'parent',           axC2, ...
        'nodeMetric',       tRates, ...
        'edgeMode',         'weight', ...
        'edgeThresholdPct', opt.edgeThresholdPct, ...
        'title',            sprintf('%s network (nodes = spikes/min)', labels.treatment));

    axC3 = nexttile(tl, 9);
    deltaRates = tRates - bRates;
    plot_network_on_mea(delta, ...
        'parent',           axC3, ...
        'nodeMetric',       deltaRates, ...
        'edgeMode',         'delta', ...
        'edgeThresholdPct', opt.edgeThresholdPct, ...
        'nodeCmap',         'turbo', ...
        'title',            '\Delta network (red = gained, blue = lost)');

    % --- Save -----------------------------------------------------------
    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end
    outFile = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_connectivity_exemplar.png', study));
    save_figure(fig, outFile);
    fprintf('fig_connectivity_exemplar(%s): saved %s\n', study, outFile);

    % --- Numeric sidecar ------------------------------------------------
    stats = struct( ...
        'study',                      study, ...
        'pair_index',                 pairIdx, ...
        'pair_count',                 numel(pairs), ...
        'baseline_dataset',           pair.baseline, ...
        'treatment_dataset',          pair.treatment, ...
        'strongest_baseline_edge_i',  iTop, ...
        'strongest_baseline_edge_j',  jTop, ...
        'strongest_baseline_weight',  bRes.adjacency(iTop, jTop), ...
        'largest_delta_edge_i',       iGain, ...
        'largest_delta_edge_j',       jGain, ...
        'largest_delta_value',        delta(iGain, jGain), ...
        'mean_baseline_weight',       mean(bRes.adjacency(triu(true(size(bRes.adjacency)),1)), 'omitnan'), ...
        'mean_treatment_weight',      mean(tRes.adjacency(triu(true(size(tRes.adjacency)),1)), 'omitnan'), ...
        'mean_delta',                 mean(delta(triu(true(size(delta)),1)), 'omitnan'), ...
        'bin_ms',                     opt.binMs, ...
        'max_lag_ms',                 opt.maxLagMs, ...
        'normalization',              opt.normalization, ...
        'edge_threshold_pct',         opt.edgeThresholdPct, ...
        'figure_file',                outFile);
    export_figure_stats(stats, fullfile(cfg.paths.figures_out, ...
        sprintf('%s_connectivity_exemplar_stats', study)));
end

% =========================================================================
function v = align_metric(vec, chUsed, channels)
% Return per-channel values aligned to `channels`; NaN for missing.
    v = nan(numel(channels), 1);
    [~, idx] = ismember(channels, chUsed(:)');
    ok = idx > 0;
    v(ok) = vec(idx(ok));
end

% =========================================================================
function [i, j] = top_edge(A)
% Return the indices of the single off-diagonal upper-triangular entry
% with the largest (signed) value. For 'delta' input use the largest
% positive entry; for baseline weights use the largest absolute peak.
    n = size(A, 1);
    mask = triu(true(n), 1);
    vals = A;
    vals(~mask) = -Inf;
    vals(isnan(vals)) = -Inf;
    [~, idx] = max(vals(:));
    [i, j] = ind2sub([n n], idx);
end

% =========================================================================
function plot_crosscorrelogram_pair(ax, bSpikes, tSpikes, ci, cj, bDur, tDur, opt, colors, labels, ~)
    binSec  = opt.binMs / 1000;

    [cB, lagsB] = pair_xcorr(bSpikes{ci}, bSpikes{cj}, binSec, bDur, opt.maxLagMs);
    [cT, lagsT] = pair_xcorr(tSpikes{ci}, tSpikes{cj}, binSec, tDur, opt.maxLagMs);

    hold(ax, 'on');
    plot(ax, lagsB, cB, 'Color', colors.baseline,  'LineWidth', 1.5);
    plot(ax, lagsT, cT, 'Color', colors.treatment, 'LineWidth', 1.5);
    xline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    xlabel(ax, 'Lag (ms)');
    ylabel(ax, 'Correlation (z)');
    xlim(ax, [-opt.maxLagMs opt.maxLagMs]);
    legend(ax, {labels.baseline, labels.treatment}, 'Location', 'best', 'Box', 'off');
    box(ax, 'on');
    hold(ax, 'off');
end

% =========================================================================
function [c, lags] = pair_xcorr(tsI, tsJ, binSec, durationSec, maxLagMs)
% Binned, z-scored cross-correlation for a single channel pair.
    nBins = max(1, floor(durationSec / binSec));
    edges = (0:nBins) * binSec;
    if isempty(tsI); xi = zeros(1, nBins); else, xi = histcounts(tsI, edges); end
    if isempty(tsJ); xj = zeros(1, nBins); else, xj = histcounts(tsJ, edges); end
    muI = mean(xi); sdI = std(xi); if sdI == 0; sdI = 1; end
    muJ = mean(xj); sdJ = std(xj); if sdJ == 0; sdJ = 1; end
    xi = (xi - muI) / sdI;
    xj = (xj - muJ) / sdJ;
    maxLagBins = round(maxLagMs / (binSec * 1000));
    [c, lagBins] = xcorr(xi, xj, maxLagBins, 'unbiased');
    c = c / nBins;
    lags = lagBins * binSec * 1000;
end

% =========================================================================
function plot_adjacency_heatmap(ax, A, ttl, symmetric)
% Square heatmap of the adjacency matrix. The colormap and colour limits
% are applied BEFORE the colorbar is created; the colorbar label is set
% via direct property assignment (cb.Label.String) to avoid the
% ylabel(cb,...) listener-callback bug that throws "Attempt to modify
% the tree during an update traversal" on some R2023b graphics trees.
    imagesc(ax, A);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');
    title(ax, ttl, 'FontWeight', 'bold');
    xlabel(ax, 'Channel index');
    ylabel(ax, 'Channel index');

    if symmetric
        finiteA = A(~isnan(A));
        absMax = max(abs(finiteA), [], 'all');
        if isempty(absMax) || absMax == 0; absMax = 1; end
        colormap(ax, diverging_cmap_local(256));
        clim(ax, [-absMax, absMax]);
    else
        finiteA = A(~isnan(A) & ~isinf(A));
        if isempty(finiteA)
            loHi = [0 1];
        else
            loHi = [min(finiteA), max(finiteA)];
            if diff(loHi) == 0; loHi = loHi + [-1 1]; end
        end
        colormap(ax, parula(256));
        clim(ax, loHi);
    end

    drawnow limitrate nocallbacks;
    cb = colorbar(ax);
    cb.Label.String     = 'Peak cross-correlation (z)';
    cb.Label.FontWeight = 'bold';
end

% =========================================================================
function cmap = diverging_cmap_local(n)
    half = floor(n / 2);
    top  = n - half;
    blue = [linspace(0.11,1,half).', linspace(0.30,1,half).', linspace(0.60,1,half).'];
    red  = [linspace(1,0.75,top).',  linspace(1,0.10,top).',  linspace(1,0.15,top).'];
    cmap = [blue; red];
end
