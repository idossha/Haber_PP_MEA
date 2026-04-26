function [results, summaryTable, stats] = fig_connectivity_summary(study, varargin)
%FIG_CONNECTIVITY_SUMMARY Main-text connectivity figure, all pairs pooled.
%
%   [results, summary, stats] = FIG_CONNECTIVITY_SUMMARY(study) runs
%   run_connectivity on every pair of STUDY ('doi' or 'ket'), extracts
%   the graph-theoretic network metrics, and renders a composite figure
%   (v1.1 layout):
%
%     Panel 1  - Paired dot plot: characteristic path length
%     Panel 2  - Paired dot plot: modularity Q
%     Panel 3  - Forest plot: effect sizes (Hedges' g) for all 6 metrics
%     Panel 4  - Spatial decay of functional connectivity
%
%   Every network-level metric is also run through paired_stats() so
%   the function returns a fully-populated stats struct array that can
%   be dropped straight into a paper table or a supplementary CSV.
%
%   FIG_CONNECTIVITY_SUMMARY(study, ...
%         'channels',      1:60,      ...
%         'binMs',          1,         ...
%         'maxLagMs',       100,       ...
%         'normalization',  'pearson', ...
%         'edgeThreshold',  0.1,       ...
%         'densitySweep',   [],        ...
%         'nullN',          50,        ...
%         'save',           false)
%   overrides defaults. If 'densitySweep' is a vector (e.g.
%   [0.05 0.10 0.15 0.20]) the function re-computes all network metrics
%   at each target density (Rubinov & Sporns 2010 recommendation) and
%   the main figure is rendered at the first value in the sweep, while
%   the whole sweep is returned in results(k).baseline.metricsSweep.
%
% INPUTS:
%   study  -  'doi' | 'ket'
%
% OUTPUTS:
%   results       -  Struct array from run_connectivity.
%   summaryTable  -  Tidy table, one row per (pair, condition).
%   stats         -  Struct array, one entry per metric, with paired_stats
%                    fields (hierarchical bootstrap CI + p, Wilcoxon,
%                    Hedges' g_av).
%
%   Also writes <study>_connectivity_summary.png to output/fig{3,5}/panels/.
%
% See also: RUN_CONNECTIVITY, NETWORK_METRICS, PAIRED_STATS,
%           FIG_CONNECTIVITY_EXEMPLAR.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',      cfg.channels.default);
    addParameter(p, 'binMs',          cfg.connectivity.bin_ms,         @(x) isscalar(x) && x > 0);
    addParameter(p, 'maxLagMs',       cfg.connectivity.max_lag_ms,     @(x) isscalar(x) && x > 0);
    addParameter(p, 'normalization',  cfg.connectivity.normalization);
    addParameter(p, 'edgeThreshold',  cfg.connectivity.edge_threshold, @(x) isempty(x) || isscalar(x));
    addParameter(p, 'edgeDensity',    cfg.connectivity.edge_density,   @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'densitySweep',   [],         @(x) isempty(x) || (isnumeric(x) && all(x > 0 & x < 1)));
    addParameter(p, 'nullN',          cfg.connectivity.null_N,         @(x) isscalar(x) && x >= 0);
    addParameter(p, 'save',           false,      @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    % --- Run connectivity for every pair (uses cached spike times) ------
    results = run_connectivity(study, ...
        'channels',      opt.channels, ...
        'binMs',         opt.binMs, ...
        'maxLagMs',      opt.maxLagMs, ...
        'normalization', opt.normalization, ...
        'edgeThreshold', opt.edgeThreshold, ...
        'edgeDensity',   opt.edgeDensity, ...
        'save',          opt.save);

    % --- Optional density sweep ----------------------------------------
    if ~isempty(opt.densitySweep)
        for k = 1:numel(results)
            results(k).baseline.metricsSweep  = compute_metrics_sweep(...
                results(k).baseline.adjacency, opt.densitySweep, opt.nullN);
            results(k).treatment.metricsSweep = compute_metrics_sweep(...
                results(k).treatment.adjacency, opt.densitySweep, opt.nullN);
        end
    end

    % --- Recompute at desired nullN if non-default ---------------------
    if opt.nullN ~= 50
        for k = 1:numel(results)
            bThresh = results(k).options.thresholdBaseline;
            tThresh = results(k).options.thresholdTreatment;
            results(k).baseline.metrics  = network_metrics(results(k).baseline.adjacency, ...
                'threshold', bThresh, 'nullN', opt.nullN);
            results(k).treatment.metrics = network_metrics(results(k).treatment.adjacency, ...
                'threshold', tThresh, 'nullN', opt.nullN);
        end
    end

    [~, labels] = get_pairs_and_labels(cfg, study);

    summaryTable = build_summary_table(results, labels);

    % All six metrics: the two featured panels + four for the forest plot.
    metricSpecs = {
        'density',              'Edge density'
        'clusteringMean',       'Clustering coeff.'
        'meanShortestPath',     'Char. path length'
        'smallWorldnessSigma',  'Small-worldness \sigma'
        'globalEfficiency',     'Global efficiency'
        'modularity',           'Modularity Q'};

    % --- paired_stats per metric ---------------------------------------
    stats = struct('metric', {}, 'label', {}, 'baselineVals', {}, ...
                   'treatmentVals', {}, 'summary', {});
    for m = 1:size(metricSpecs, 1)
        field = metricSpecs{m, 1};
        bVals = arrayfun(@(r) safe_get(r.baseline.metrics,  field), results);
        tVals = arrayfun(@(r) safe_get(r.treatment.metrics, field), results);
        stats(m).metric        = field;
        stats(m).label         = metricSpecs{m, 2};
        stats(m).baselineVals  = bVals(:);
        stats(m).treatmentVals = tVals(:);
        stats(m).summary = paired_stats(bVals, tVals);
    end

    % BH-FDR across the family of network metrics (reviewer courtesy).
    pFamily = arrayfun(@(s) s.summary.bootstrap.pDelta, stats);
    [pAdj, rejected] = bh_fdr(pFamily, 0.05);
    for m = 1:numel(stats)
        stats(m).summary.bootstrap.pDeltaBHFDR = pAdj(m);
        stats(m).summary.bootstrap.significantBH = rejected(m);
    end

    % --- Figure layout (v1.1: 1x4 combined) ----------------------------
    % Panel 1: path length paired | Panel 2: modularity paired
    % Panel 3: forest plot (all 6 metrics) | Panel 4: spatial decay
    colors = paired_plot_colors(study);
    fig = create_panel_figure(18.0, 3.6);
    tl = tiledlayout(fig, 1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

    idxPL  = find(strcmp({stats.metric}, 'meanShortestPath'));
    idxMod = find(strcmp({stats.metric}, 'modularity'));

    ax1 = nexttile(tl, 1);
    plot_paired_panel(ax1, stats(idxPL), labels, colors);

    ax2 = nexttile(tl, 2);
    plot_paired_panel(ax2, stats(idxMod), labels, colors);

    ax3 = nexttile(tl, 3);
    plot_forest(ax3, stats, colors);

    ax4 = nexttile(tl, 4);
    plot_spatial_decay(ax4, results, colors, labels);

    % --- Spatial-decay slope test (per-well regression slope) -------------
    slopeStats = spatial_decay_slope_test(results);
    if slopeStats.pSlope < 0.001
        pStr = 'p_{slope} < 0.001';
    else
        pStr = sprintf('p_{slope} = %.3f', slopeStats.pSlope);
    end
    text(ax4, 0.97, 0.97, sprintf('Slope test: %s', pStr), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 7, 'FontName', 'Arial');

    % --- Nature/NPP styling -----------------------------------------------
    apply_nature_style(fig);

    panelDir = output_path(cfg, study, 'connectivity', 'panels');
    statsDir = output_path(cfg, study, 'connectivity', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    outFile = fullfile(panelDir, ...
        sprintf('%s_connectivity_summary.png', study));
    save_figure(fig, outFile);
    fprintf('fig_connectivity_summary(%s): saved %s\n', study, outFile);

    % --- Numeric sidecar: one row per metric ---------------------------
    metricTable = struct('metric', {}, 'label', {}, ...
        'median_baseline', {}, 'median_treatment', {}, ...
        'median_delta', {}, ...
        'ci_delta_lo', {}, 'ci_delta_hi', {}, ...
        'p_bootstrap', {}, 'p_bh_fdr', {}, ...
        'p_wilcoxon', {}, 'hedges_g_av', {});
    for m = 1:numel(stats)
        s = stats(m).summary;
        row = struct();
        row.metric            = stats(m).metric;
        row.label             = stats(m).label;
        row.median_baseline   = s.medianBaseline;
        row.median_treatment  = s.medianTreatment;
        row.median_delta      = s.medianDelta;
        row.ci_delta_lo       = s.bootstrap.ciDelta(1);
        row.ci_delta_hi       = s.bootstrap.ciDelta(2);
        row.p_bootstrap       = s.bootstrap.pDelta;
        if isfield(s.bootstrap, 'pDeltaBHFDR')
            row.p_bh_fdr      = s.bootstrap.pDeltaBHFDR;
        else
            row.p_bh_fdr      = NaN;
        end
        row.p_wilcoxon        = s.wilcoxon.p;
        row.hedges_g_av       = s.hedgesGav;
        metricTable(end+1) = row; %#ok<AGROW>
    end
    % One stats file per metric (CSV/JSON).
    for m = 1:numel(metricTable)
        export_figure_stats(metricTable(m), ...
            fullfile(statsDir, sprintf('%s_%s', study, metricTable(m).metric)));
    end
    % Tidy-table CSV for paper-writing.
    tidyPath = fullfile(statsDir, ...
        sprintf('%s_connectivity_metrics.csv', study));
    try
        writetable(struct2table(metricTable), tidyPath);
    catch ME
        warning('fig_connectivity_summary:CsvFailed', ...
            'Could not write tidy table: %s', ME.message);
    end
    fprintf('fig_connectivity_summary(%s): stats at %s\n', study, statsDir);
end

% =========================================================================
function v = safe_get(s, field)
    if isfield(s, field); v = s.(field); else, v = NaN; end
end

% =========================================================================
function metricsSweep = compute_metrics_sweep(adjacency, densities, nullN)
% Run network_metrics at a vector of target densities by taking the
% appropriate percentile threshold of the weight distribution.
    n = size(adjacency, 1);
    upperMask = triu(true(n), 1);
    vals = adjacency(upperMask);
    vals = vals(~isnan(vals));
    metricsSweep = struct('density', {}, 'threshold', {}, 'metrics', {});
    for d = 1:numel(densities)
        targetDensity = densities(d);
        if isempty(vals)
            thresh = 0;
        else
            pct = 100 * (1 - targetDensity);
            thresh = prctile(vals, pct);
        end
        metricsSweep(d).density   = targetDensity;
        metricsSweep(d).threshold = thresh;
        metricsSweep(d).metrics   = network_metrics(adjacency, ...
            'threshold', thresh, 'nullN', nullN);
    end
end

% =========================================================================
function T = build_summary_table(results, labels)
    nPairs = numel(results);
    rows = cell(2 * nPairs, 1);
    fields = {'density','degreeMean','clusteringMean','meanShortestPath', ...
              'globalEfficiency','modularity','nCommunities','smallWorldnessSigma'};
    idx = 1;
    for k = 1:nPairs
        for side = {'baseline','treatment'}
            s = results(k).(side{1});
            m = s.metrics;
            row = struct();
            row.pair        = k;
            row.condition   = side{1};
            row.label       = labels.(side{1});
            row.datasetName = s.datasetName;
            for f = 1:numel(fields)
                if isfield(m, fields{f})
                    row.(fields{f}) = m.(fields{f});
                else
                    row.(fields{f}) = NaN;
                end
            end
            rows{idx} = row;
            idx = idx + 1;
        end
    end
    T = struct2table(vertcat(rows{:}));
end

% =========================================================================
function plot_paired_panel(ax, stat, labels, colors)
%PLOT_PAIRED_PANEL  Refined paired dot plot with IQR shading and bracket.
    hold(ax, 'on');
    bVals = stat.baselineVals;
    tVals = stat.treatmentVals;
    n = numel(bVals);

    % Deterministic jitter based on rank order
    jit = linspace(-0.04, 0.04, n)';

    % IQR shading
    bIQR = [prctile(bVals, 25), prctile(bVals, 75)];
    tIQR = [prctile(tVals, 25), prctile(tVals, 75)];
    fill(ax, [0.7 1.3 1.3 0.7], [bIQR(1) bIQR(1) bIQR(2) bIQR(2)], ...
        colors.baseline, 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    fill(ax, [1.7 2.3 2.3 1.7], [tIQR(1) tIQR(1) tIQR(2) tIQR(2)], ...
        colors.treatment, 'FaceAlpha', 0.12, 'EdgeColor', 'none');

    % Paired connecting lines
    for k = 1:n
        if isnan(bVals(k)) || isnan(tVals(k)); continue; end
        if tVals(k) > bVals(k)
            lc = [colors.increase 0.6];
        else
            lc = [colors.decrease 0.6];
        end
        plot(ax, [1+jit(k) 2+jit(k)], [bVals(k) tVals(k)], '-', ...
            'Color', lc, 'LineWidth', 1.4);
    end

    % Dots
    scatter(ax, 1+jit, bVals, 70, colors.baseline, 'filled', ...
        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.6);
    scatter(ax, 2+jit, tVals, 70, colors.treatment, 'filled', ...
        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.6);

    % Median diamonds
    plot(ax, 1, median(bVals, 'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(ax, 2, median(tVals, 'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

    % Significance bracket
    yMax = max([bVals; tVals], [], 'omitnan');
    yBracket = yMax * 1.08;
    sigBH = isfield(stat.summary.bootstrap, 'significantBH') && ...
            stat.summary.bootstrap.significantBH;
    if sigBH; star = '*'; else; star = 'n.s.'; end
    plot(ax, [1 2], [yBracket yBracket], 'k-', 'LineWidth', 0.8);
    text(ax, 1.5, yBracket * 1.01, star, ...
        'HorizontalAlignment', 'center', 'FontSize', 13, ...
        'FontName', 'Arial', 'FontWeight', 'bold');

    % P-value and effect size below bracket
    pBoot = stat.summary.bootstrap.pDelta;
    if pBoot < 0.001; pStr = 'p < .001';
    else; pStr = sprintf('p = %.3f', pBoot); end
    gStr = sprintf('g = %.2f', stat.summary.hedgesGav);
    text(ax, 1.5, yBracket * 1.07, {pStr, gStr}, ...
        'HorizontalAlignment', 'center', 'FontSize', 8, ...
        'FontName', 'Arial', 'Color', [0.35 0.35 0.35]);

    set(ax, 'XTick', [1 2], 'XTickLabel', {labels.baseline, labels.treatment});
    xlim(ax, [0.4 2.6]);
    yPad = (yBracket * 1.12 - min([bVals; tVals], [], 'omitnan')) * 0.02;
    ylim(ax, [min([bVals; tVals], [], 'omitnan') - yPad, yBracket * 1.12]);
    ylabel(ax, stat.label, 'FontSize', 11, 'FontWeight', 'bold');
    set(ax, 'FontSize', 10, 'FontName', 'Arial', 'TickDir', 'out', ...
        'LineWidth', 1, 'Box', 'off');
    hold(ax, 'off');
end

% =========================================================================
function plot_forest(ax, statsArr, colors)
%PLOT_FOREST  Horizontal forest plot of effect sizes for all metrics.
    hold(ax, 'on');
    nMetrics = numel(statsArr);

    gVals    = arrayfun(@(s) s.summary.hedgesGav, statsArr);
    medDelta = arrayfun(@(s) s.summary.medianDelta, statsArr);
    ciLo     = arrayfun(@(s) s.summary.bootstrap.ciDelta(1), statsArr);
    ciHi     = arrayfun(@(s) s.summary.bootstrap.ciDelta(2), statsArr);

    % Convert delta CI to g-scale
    ciLoG = nan(1, nMetrics);
    ciHiG = nan(1, nMetrics);
    for m = 1:nMetrics
        if abs(gVals(m)) > 0.001 && abs(medDelta(m)) > eps
            pooledSD = abs(medDelta(m)) / abs(gVals(m));
        else
            pooledSD = std([statsArr(m).baselineVals; statsArr(m).treatmentVals], 'omitnan');
            if pooledSD < eps; pooledSD = 1; end
        end
        ciLoG(m) = ciLo(m) / pooledSD;
        ciHiG(m) = ciHi(m) / pooledSD;
    end

    yPos = 1:nMetrics;

    % Zero reference line
    xline(ax, 0, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);

    % Negligible-effect zone (|g| < 0.2)
    fill(ax, [-0.2 0.2 0.2 -0.2], [0 0 nMetrics+1 nMetrics+1], ...
        [0.96 0.96 0.96], 'EdgeColor', 'none');

    for m = 1:nMetrics
        y = yPos(m);
        g = gVals(m);
        sigBH = isfield(statsArr(m).summary.bootstrap, 'significantBH') && ...
                statsArr(m).summary.bootstrap.significantBH;

        % CI whisker
        plot(ax, [ciLoG(m) ciHiG(m)], [y y], '-', ...
            'Color', [0.45 0.45 0.45], 'LineWidth', 1.8);

        % Marker: diamond if significant, circle if not
        if sigBH
            scatter(ax, g, y, 110, colors.treatment, 'filled', 'd', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1);
        else
            scatter(ax, g, y, 70, [0.6 0.6 0.6], 'filled', 'o', ...
                'MarkerEdgeColor', [0.4 0.4 0.4]);
        end
    end

    set(ax, 'YTick', yPos, 'YTickLabel', {statsArr.label});
    ylim(ax, [0.4 nMetrics + 0.6]);
    xlabel(ax, 'Hedges'' g_{av}', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax, 'Effect sizes (all metrics)', 'FontSize', 11, 'FontWeight', 'bold');
    set(ax, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', ...
        'LineWidth', 1, 'Box', 'off');

    % Legend
    h1 = scatter(ax, NaN, NaN, 90, colors.treatment, 'filled', 'd', 'MarkerEdgeColor', 'k');
    h2 = scatter(ax, NaN, NaN, 60, [0.6 0.6 0.6], 'filled', 'o', 'MarkerEdgeColor', [0.4 0.4 0.4]);
    legend(ax, [h1 h2], {'BH-FDR sig.', 'Not sig.'}, ...
        'Location', 'southeast', 'Box', 'off', 'FontSize', 8);

    hold(ax, 'off');
end

% =========================================================================
function plot_spatial_decay(ax, results, colors, labels)
% Binned-median of |weight| vs inter-electrode distance, pooled across pairs.
    layout = mea60_layout();
    D = layout.distanceMatrix;
    upperMask = triu(true(size(D)), 1);
    distVec = D(upperMask);

    bW = []; tW = [];
    for k = 1:numel(results)
        bA = results(k).baseline.adjacency;
        tA = results(k).treatment.adjacency;
        bA(isnan(bA)) = 0;
        tA(isnan(tA)) = 0;
        bW = [bW; bA(upperMask)]; %#ok<AGROW>
        tW = [tW; tA(upperMask)]; %#ok<AGROW>
    end

    % Match distance vector to pooled weights.
    repDist = repmat(distVec, numel(results), 1);

    hold(ax, 'on');
    edges = 0:1:8;
    [bMed, bLo, bHi] = binned_ci(repDist, abs(bW), edges);
    [tMed, tLo, tHi] = binned_ci(repDist, abs(tW), edges);
    ctrs = edges(1:end-1) + diff(edges)/2;

    fill(ax, [ctrs fliplr(ctrs)], [bLo fliplr(bHi)], colors.baseline, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    fill(ax, [ctrs fliplr(ctrs)], [tLo fliplr(tHi)], colors.treatment, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax, ctrs, bMed, 'o-', 'Color', colors.baseline,  'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', colors.baseline);
    plot(ax, ctrs, tMed, 's-', 'Color', colors.treatment, 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', colors.treatment);

    xlabel(ax, 'Inter-electrode distance (100 \mum)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax, '|Cross-correlation peak|', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax, 'Spatial decay', 'FontSize', 11, 'FontWeight', 'bold');
    legend(ax, {labels.baseline, labels.treatment}, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    set(ax, 'FontSize', 9, 'FontName', 'Arial', 'TickDir', 'out', ...
        'LineWidth', 1, 'Box', 'off');
    hold(ax, 'off');
end

% =========================================================================
function [m, lo, hi] = binned_ci(x, y, edges)
    m = nan(1, numel(edges)-1);
    lo = m; hi = m;
    for b = 1:numel(edges)-1
        in = x >= edges(b) & x < edges(b+1);
        v = y(in);
        v = v(~isnan(v));
        if isempty(v); continue; end
        m(b)  = median(v);
        lo(b) = prctile(v, 25);
        hi(b) = prctile(v, 75);
    end
end

% =========================================================================
function s = format_p(p)
    if isnan(p); s = '  n/a '; return; end
    if p < 0.001
        s = '<.001';
    else
        s = sprintf('%.3f', p);
    end
end

% =========================================================================
function out = spatial_decay_slope_test(results)
%SPATIAL_DECAY_SLOPE_TEST  Per-well OLS slope of |weight| ~ distance.
%   Tests whether the spatial decay slope differs between baseline and
%   treatment using a paired Wilcoxon signed-rank test on per-well slopes.
    layout   = mea60_layout();
    D        = layout.distanceMatrix;
    upper    = triu(true(size(D)), 1);
    distVec  = D(upper);
    nPairs   = numel(results);

    baseSlopes  = nan(nPairs, 1);
    treatSlopes = nan(nPairs, 1);
    for k = 1:nPairs
        bA = results(k).baseline.adjacency;
        tA = results(k).treatment.adjacency;
        bA(isnan(bA)) = 0;
        tA(isnan(tA)) = 0;
        baseSlopes(k)  = ols_slope(abs(bA(upper)), distVec);
        treatSlopes(k) = ols_slope(abs(tA(upper)), distVec);
    end

    slopeDiff  = treatSlopes - baseSlopes;
    if nPairs >= 5
        [pSlope, ~, w] = signrank(baseSlopes, treatSlopes);
    else
        [pSlope, ~, w] = signrank(baseSlopes, treatSlopes, 'method', 'exact');
    end

    out.baseSlopes  = baseSlopes;
    out.treatSlopes = treatSlopes;
    out.slopeDiff   = slopeDiff;
    out.medianDiff  = median(slopeDiff);
    out.pSlope      = pSlope;
    if isfield(w, 'signedrank')
        out.statistic = w.signedrank;
    else
        out.statistic = NaN;
    end
end

% =========================================================================
function slope = ols_slope(w, d)
%OLS_SLOPE  Simple OLS slope of w on d.
    X = [ones(numel(d), 1), d(:)];
    b = X \ w(:);
    slope = b(2);
end
