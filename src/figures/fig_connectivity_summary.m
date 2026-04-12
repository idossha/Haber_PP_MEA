function [results, summaryTable, stats] = fig_connectivity_summary(study, varargin)
%FIG_CONNECTIVITY_SUMMARY Main-text connectivity figure, all pairs pooled.
%
%   [results, summary, stats] = FIG_CONNECTIVITY_SUMMARY(study) runs
%   run_connectivity on every pair of STUDY ('doi' or 'ket'), extracts
%   the graph-theoretic network metrics, and renders a composite figure
%   following the Sharf 2022 / Downes 2012 / Varley 2024 template:
%
%     Row 1  - Network-metric slope plots (one panel per metric,
%              paired per-well lines + group median + 95% bootstrap CI)
%     Row 2  - Spatial decay  |  edge-weight KDE
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
%   Also writes <study>_connectivity_summary.png to cfg.paths.figures_out.
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

    metricSpecs = {
        'density',              'Edge density'
        'clusteringMean',       'Clustering coefficient'
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

    % --- Figure layout --------------------------------------------------
    colors = paired_plot_colors(study);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 1000]);
    tl = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Connectivity effects: %s vs %s (n = %d pairs)', ...
          labels.baseline, labels.treatment, numel(results)));

    % Row 1-2 : 6 network-metric slope plots
    for m = 1:numel(stats)
        ax = nexttile(tl, m);
        plot_paired_metric(ax, stats(m), colors, labels);
    end

    % Row 3 : spatial decay + edge-weight histogram + (empty for padding)
    axDecay = nexttile(tl, 7);
    plot_spatial_decay(axDecay, results, colors, labels);

    axHist = nexttile(tl, 8);
    plot_edge_weight_kde(axHist, results, colors, labels);

    axLegend = nexttile(tl, 9);
    render_summary_text(axLegend, stats, opt);

    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end
    outFile = fullfile(cfg.paths.figures_out, ...
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
    % One stats file per metric (CSV/JSON) and a single summary table.
    baseDir = fullfile(cfg.paths.figures_out, sprintf('%s_connectivity_stats', study));
    if ~exist(baseDir, 'dir'); mkdir(baseDir); end
    for m = 1:numel(metricTable)
        export_figure_stats(metricTable(m), ...
            fullfile(baseDir, sprintf('%s_%s', study, metricTable(m).metric)));
    end
    % Tidy-table CSV for paper-writing.
    tidyPath = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_connectivity_metrics.csv', study));
    try
        writetable(struct2table(metricTable), tidyPath);
    catch ME
        warning('fig_connectivity_summary:CsvFailed', ...
            'Could not write tidy table: %s', ME.message);
    end
    fprintf('fig_connectivity_summary(%s): stats at %s\n', study, baseDir);
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
function plot_paired_metric(ax, stat, colors, labels)
    bVals = stat.baselineVals;
    tVals = stat.treatmentVals;
    n = numel(bVals);

    hold(ax, 'on');
    for k = 1:n
        if isnan(bVals(k)) || isnan(tVals(k)); continue; end
        if tVals(k) > bVals(k); col = [0.1 0.55 0.25];
        elseif tVals(k) < bVals(k); col = [0.75 0.20 0.20];
        else; col = [0.7 0.7 0.7]; end
        plot(ax, [1 2], [bVals(k) tVals(k)], '-', ...
            'Color', [col 0.6], 'LineWidth', 1);
    end
    scatter(ax, ones(n,1),    bVals, 60, 'filled', ...
        'MarkerFaceColor', colors.baseline,  'MarkerEdgeColor', 'k');
    scatter(ax, 2*ones(n,1),  tVals, 60, 'filled', ...
        'MarkerFaceColor', colors.treatment, 'MarkerEdgeColor', 'k');

    if any(~isnan(bVals))
        plot(ax, 1, mean(bVals, 'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    end
    if any(~isnan(tVals))
        plot(ax, 2, mean(tVals, 'omitnan'), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    end

    % Annotate significance marker from BH-adjusted hierarchical bootstrap.
    yMax = max([bVals; tVals], [], 'omitnan');
    if ~isempty(stat.summary) && ~isnan(stat.summary.bootstrap.pDelta)
        pStr = format_p(stat.summary.bootstrap.pDelta);
        if isfield(stat.summary.bootstrap, 'pDeltaBHFDR') && stat.summary.bootstrap.significantBH
            pStr = ['*' pStr];
        end
        text(ax, 1.5, yMax, pStr, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    end

    set(ax, 'XTick', [1 2], 'XTickLabel', {labels.baseline, labels.treatment});
    xlim(ax, [0.5 2.5]);
    ylabel(ax, stat.label);
    box(ax, 'on');
    hold(ax, 'off');
end

% =========================================================================
function plot_spatial_decay(ax, results, colors, labels)
% Scatter and binned-median of |weight| vs inter-electrode distance,
% pooled across pairs.
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
    [bMean, bLo, bHi] = binned_ci(repDist, abs(bW), edges);
    [tMean, tLo, tHi] = binned_ci(repDist, abs(tW), edges);
    ctrs = edges(1:end-1) + diff(edges)/2;

    fill(ax, [ctrs fliplr(ctrs)], [bLo fliplr(bHi)], colors.baseline, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    fill(ax, [ctrs fliplr(ctrs)], [tLo fliplr(tHi)], colors.treatment, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax, ctrs, bMean, '-', 'Color', colors.baseline,  'LineWidth', 2);
    plot(ax, ctrs, tMean, '-', 'Color', colors.treatment, 'LineWidth', 2);

    xlabel(ax, 'Inter-electrode distance (100 \mum units)');
    ylabel(ax, '|peak cross-correlation|  (log scale)');
    title(ax, 'Spatial decay of functional connectivity', 'FontWeight', 'bold');
    legend(ax, {[labels.baseline ' IQR'], [labels.treatment ' IQR'], ...
                labels.baseline, labels.treatment}, ...
        'Location', 'northeast', 'Box', 'off');
    % Log Y reveals the exponential decay otherwise hidden by the near-zero
    % floor.
    set(ax, 'YScale', 'log');
    posMask = [bMean tMean bLo tLo bHi tHi];
    posMask = posMask(posMask > 0);
    if ~isempty(posMask)
        lo = min(posMask) * 0.7;
        hi = max(posMask) * 1.3;
        ylim(ax, [lo, hi]);
    end
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
function plot_edge_weight_kde(ax, results, colors, labels)
% Kernel-density estimates of the pooled edge-weight distribution for
% baseline vs treatment.
    bW = []; tW = [];
    for k = 1:numel(results)
        bA = results(k).baseline.adjacency;
        tA = results(k).treatment.adjacency;
        upperMask = triu(true(size(bA)), 1);
        bW = [bW; bA(upperMask)]; %#ok<AGROW>
        tW = [tW; tA(upperMask)]; %#ok<AGROW>
    end
    bW = bW(~isnan(bW));
    tW = tW(~isnan(tW));

    if isempty(bW) || isempty(tW)
        text(ax, 0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return;
    end

    % Log-X axis: negative and zero edges are dropped for the KDE.
    bPos = bW(bW > 0);
    tPos = tW(tW > 0);
    if isempty(bPos) || isempty(tPos)
        text(ax, 0.5, 0.5, 'No positive edges', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return;
    end
    xMin = min([bPos; tPos]);
    xMax = max([bPos; tPos]);
    % KDE in log space: evaluate density at log10(weight).
    pts = logspace(log10(xMin), log10(xMax), 200);
    fB = ksdensity(log10(bPos), log10(pts));
    fT = ksdensity(log10(tPos), log10(pts));

    hold(ax, 'on');
    fill(ax, [pts fliplr(pts)], [fB zeros(size(pts))], colors.baseline, ...
        'FaceAlpha', 0.35, 'EdgeColor', colors.baseline, 'LineWidth', 1.5);
    fill(ax, [pts fliplr(pts)], [fT zeros(size(pts))], colors.treatment, ...
        'FaceAlpha', 0.35, 'EdgeColor', colors.treatment, 'LineWidth', 1.5);
    xlabel(ax, 'Edge weight  (log scale)');
    ylabel(ax, 'Density (log_{10} weight)');
    title(ax, 'Edge-weight distribution', 'FontWeight', 'bold');
    set(ax, 'XScale', 'log');
    legend(ax, {labels.baseline, labels.treatment}, ...
        'Location', 'northeast', 'Box', 'off');
    hold(ax, 'off');
end

% =========================================================================
function render_summary_text(ax, stats, opt)
    axis(ax, 'off');
    lines = {};
    lines{end+1} = sprintf('\\bf Statistics summary');
    lines{end+1} = '';
    lines{end+1} = sprintf('Bin = %g ms', opt.binMs);
    lines{end+1} = sprintf('Max lag = %g ms', opt.maxLagMs);
    if isempty(opt.edgeThreshold)
        lines{end+1} = sprintf('Edge density target = %.2f', opt.edgeDensity);
    else
        lines{end+1} = sprintf('Edge threshold (abs.) = %.2g', opt.edgeThreshold);
    end
    lines{end+1} = '';
    lines{end+1} = 'Metric  |  pBoot  |  pBH  |  g_{av}';
    for m = 1:numel(stats)
        s = stats(m).summary;
        pB = s.bootstrap.pDelta;
        pA = s.bootstrap.pDeltaBHFDR;
        g  = s.hedgesGav;
        lines{end+1} = sprintf('%s  |  %s  |  %s  |  %.2f', ...
            stats(m).label, format_p(pB), format_p(pA), g); %#ok<AGROW>
    end
    y = 0.97;
    dy = 0.065;
    for k = 1:numel(lines)
        text(ax, 0.02, y, lines{k}, 'FontName', 'Menlo', 'FontSize', 9, ...
            'Interpreter', 'tex');
        y = y - dy;
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
