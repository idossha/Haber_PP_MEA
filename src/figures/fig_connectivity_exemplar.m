function fig_connectivity_exemplar(study, varargin)
%FIG_CONNECTIVITY_EXEMPLAR Individual connectivity panels for one pair.
%
%   FIG_CONNECTIVITY_EXEMPLAR(study) generates individual panel files for
%   one representative pair from STUDY ('doi' or 'ket'):
%
%     2 correlogram CSVs   (strongest baseline edge, largest delta edge)
%                          — plotted by scripts/plot_correlograms.py
%     3 adjacency heatmaps (baseline, treatment, delta) — spatially ordered
%     6 network graphs     (baseline, treatment, delta) x (clean, labeled)
%
%   Each panel is saved as its own 600 DPI PDF + PNG. Correlogram data is
%   exported as CSV + JSON metadata for Python rendering.
%
%   Output subdirectories under the figure root:
%     correlogram/   heatmap/   network/
%
% See also: CONNECTIVITY_XCORR, PLOT_NETWORK_ON_MEA, RUN_FIGURES.

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
    addParameter(p, 'outDir',          output_path(cfg, study, 'connectivity', ''));
    addParameter(p, 'statsDir',        output_path(cfg, study, 'connectivity', 'stats'));
    addParameter(p, 'prefix',          '');
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    % --- Pair selection ---------------------------------------------------
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

    if isempty(opt.prefix)
        prefix = study;
    else
        prefix = opt.prefix;
    end

    fprintf('fig_connectivity_exemplar(%s): pair %d/%d  [prefix=%s]\n', ...
        study, pairIdx, numel(pairs), prefix);

    % --- Load spikes and rates --------------------------------------------
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, opt.channels, cfg);
    bCache = load_cache(pair.baseline,  cfg);
    tCache = load_cache(pair.treatment, cfg);

    xcorrArgs = {'binMs', opt.binMs, 'maxLagMs', opt.maxLagMs, ...
                 'normalization', opt.normalization};
    bRes = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, xcorrArgs{:});
    tRes = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, xcorrArgs{:});
    delta = tRes.adjacency - bRes.adjacency;

    bRates     = align_metric(bCache.spikeRates, bCache.channelsUsed, opt.channels);
    tRates     = align_metric(tCache.spikeRates, tCache.channelsUsed, opt.channels);
    deltaRates = tRates - bRates;

    % --- Output subdirectories --------------------------------------------
    % Correlograms go to SI (supplementary) — they are quality-control
    % panels, not a main-figure result (Cutts & Eglen 2014 §1.6).
    siDir   = output_path(cfg, '', 'si', '');
    corrDir = fullfile(siDir, 'correlogram');
    heatDir = fullfile(opt.outDir, 'heatmap');
    netDir  = fullfile(opt.outDir, 'network');
    for d = {corrDir, heatDir, netDir, opt.statsDir}
        if ~exist(d{1}, 'dir'); mkdir(d{1}); end
    end

    colors = paired_plot_colors(study);
    layout = mea60_layout();

    % --- Exemplar edges ---------------------------------------------------
    [iTop, jTop]   = top_edge(bRes.adjacency);
    [iGain, jGain] = top_edge(delta);

    % === Correlograms (export CSV + JSON for Python plotting) =============
    export_correlogram(bSpikes, tSpikes, iTop, jTop, bMeta, tMeta, opt, ...
        labels, colors, layout, ...
        sprintf('Strongest baseline edge (el. %d - %d)', ...
            layout.mcsLabels(iTop), layout.mcsLabels(jTop)), ...
        fullfile(corrDir, [prefix '_correlogram_top_baseline']));

    export_correlogram(bSpikes, tSpikes, iGain, jGain, bMeta, tMeta, opt, ...
        labels, colors, layout, ...
        sprintf('Largest delta edge (el. %d - %d)', ...
            layout.mcsLabels(iGain), layout.mcsLabels(jGain)), ...
        fullfile(corrDir, [prefix '_correlogram_top_delta']));

    % === Adjacency heatmaps (spatially ordered by MCS label) ==============
    % Compute dead-channel mask once from baseline, apply to all three.
    refIdx   = cfg.channels.reference;   % PZ5 ch 15 = internal reference
    deadMask = detect_dead_channels(bRes.adjacency, layout, refIdx);

    save_adj_heatmap(bRes.adjacency, labels.baseline, false, layout, deadMask, ...
        fullfile(heatDir, [prefix '_adj_baseline']));

    save_adj_heatmap(tRes.adjacency, labels.treatment, false, layout, deadMask, ...
        fullfile(heatDir, [prefix '_adj_treatment']));

    save_adj_heatmap(delta, '\Delta (treatment - baseline)', true, layout, deadMask, ...
        fullfile(heatDir, [prefix '_adj_delta']));

    % === Network graphs (clean + labeled versions) ========================
    for mode = {'clean', 'labeled'}
        showLbl = strcmp(mode{1}, 'labeled');
        if showLbl; suffix = '_labeled'; else; suffix = ''; end

        save_network_panel(bRes.adjacency, bRates, opt, ...
            sprintf('%s network', labels.baseline), 'weight', 'parula', showLbl, ...
            fullfile(netDir, [prefix '_network_baseline' suffix]));

        save_network_panel(tRes.adjacency, tRates, opt, ...
            sprintf('%s network', labels.treatment), 'weight', 'parula', showLbl, ...
            fullfile(netDir, [prefix '_network_treatment' suffix]));

        save_network_panel(delta, deltaRates, opt, ...
            '\Delta network (red = gained, blue = lost)', 'delta', 'turbo', showLbl, ...
            fullfile(netDir, [prefix '_network_delta' suffix]));
    end

    % --- Numeric sidecar --------------------------------------------------
    sidecar = struct( ...
        'study',                      study, ...
        'pair_index',                 pairIdx, ...
        'pair_count',                 numel(pairs), ...
        'baseline_dataset',           pair.baseline, ...
        'treatment_dataset',          pair.treatment, ...
        'strongest_baseline_edge_i',  layout.mcsLabels(iTop), ...
        'strongest_baseline_edge_j',  layout.mcsLabels(jTop), ...
        'strongest_baseline_weight',  bRes.adjacency(iTop, jTop), ...
        'largest_delta_edge_i',       layout.mcsLabels(iGain), ...
        'largest_delta_edge_j',       layout.mcsLabels(jGain), ...
        'largest_delta_value',        delta(iGain, jGain), ...
        'mean_baseline_weight',       mean(bRes.adjacency(triu(true(size(bRes.adjacency)),1)), 'omitnan'), ...
        'mean_treatment_weight',      mean(tRes.adjacency(triu(true(size(tRes.adjacency)),1)), 'omitnan'), ...
        'mean_delta',                 mean(delta(triu(true(size(delta)),1)), 'omitnan'), ...
        'bin_ms',                     opt.binMs, ...
        'max_lag_ms',                 opt.maxLagMs, ...
        'normalization',              opt.normalization, ...
        'edge_threshold_pct',         opt.edgeThresholdPct);
    export_figure_stats(sidecar, ...
        fullfile(opt.statsDir, [prefix '_connectivity_exemplar_stats']));

    fprintf('  -> heatmap+network saved to %s/\n', opt.outDir);
    fprintf('  -> correlograms  saved to %s/\n', corrDir);
end


% =========================================================================
%                        HELPER FUNCTIONS
% =========================================================================

function v = align_metric(vec, chUsed, channels)
    v = nan(numel(channels), 1);
    [~, idx] = ismember(channels, chUsed(:)');
    ok = idx > 0;
    v(ok) = vec(idx(ok));
end

function [i, j] = top_edge(A)
    n = size(A, 1);
    mask = triu(true(n), 1);
    vals = A;
    vals(~mask) = -Inf;
    vals(isnan(vals)) = -Inf;
    [~, idx] = max(vals(:));
    [i, j] = ind2sub([n n], idx);
end

% -------------------------------------------------------------------------
function export_correlogram(bSpikes, tSpikes, ci, cj, bMeta, tMeta, ...
                            opt, labels, colors, layout, ttl, outBase)
%EXPORT_CORRELOGRAM Write correlogram data as CSV + JSON for Python.
    binSec = opt.binMs / 1000;
    [cB, lagsB] = pair_xcorr(bSpikes{ci}, bSpikes{cj}, binSec, ...
                             bMeta.durationSec, opt.maxLagMs);
    [cT, lagsT] = pair_xcorr(tSpikes{ci}, tSpikes{cj}, binSec, ...
                             tMeta.durationSec, opt.maxLagMs);

    % CSV: lag_ms, baseline, treatment
    T = table(lagsB(:), cB(:), cT(:), ...
        'VariableNames', {'lag_ms', 'baseline', 'treatment'});
    writetable(T, [outBase '_data.csv']);

    % JSON metadata
    meta = struct( ...
        'title',           ttl, ...
        'baseline_label',  labels.baseline, ...
        'treatment_label', labels.treatment, ...
        'baseline_color',  colors.baseline, ...
        'treatment_color', colors.treatment, ...
        'channel_i_mcs',   layout.mcsLabels(ci), ...
        'channel_j_mcs',   layout.mcsLabels(cj), ...
        'max_lag_ms',      opt.maxLagMs);
    fid = fopen([outBase '_meta.json'], 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
end

function [c, lags] = pair_xcorr(tsI, tsJ, binSec, durationSec, maxLagMs)
    nBins = max(1, floor(durationSec / binSec));
    edges = (0:nBins) * binSec;
    if isempty(tsI); xi = zeros(1, nBins); else; xi = histcounts(tsI, edges); end
    if isempty(tsJ); xj = zeros(1, nBins); else; xj = histcounts(tsJ, edges); end
    muI = mean(xi); sdI = std(xi); if sdI == 0; sdI = 1; end
    muJ = mean(xj); sdJ = std(xj); if sdJ == 0; sdJ = 1; end
    xi = (xi - muI) / sdI;
    xj = (xj - muJ) / sdJ;
    maxLagBins = round(maxLagMs / (binSec * 1000));
    [c, lagBins] = xcorr(xi, xj, maxLagBins, 'unbiased');
    lags = lagBins * binSec * 1000;
end

% -------------------------------------------------------------------------
function deadMask = detect_dead_channels(A, layout, refIdx)
%DETECT_DEAD_CHANNELS  Identify reference + dead amplifier channels.
%   Returns a logical mask (1 = exclude) over the UNSORTED channel indices.
%   Dead channels are those whose mean cross-correlation is more than 3 SD
%   below the channel mean — an extreme threshold that catches only truly
%   non-contributing hardware faults.
    nCh = size(A, 1);
    deadMask = false(1, nCh);
    deadMask(refIdx) = true;

    chMean = nan(1, nCh);
    for c = 1:nCh
        row = A(c, [1:c-1, c+1:nCh]);
        chMean(c) = mean(row, 'omitnan');
    end
    % Channels with mean xcorr below 1/10th of the median are dead
    % hardware (amplifier faults, broken connector pins). This is an
    % extreme threshold: only catches channels ~10x weaker than typical.
    med = median(chMean, 'omitnan');
    deadMask(chMean < med / 10) = true;

    excluded = layout.mcsLabels(deadMask);
    if any(deadMask)
        fprintf('    heatmap: excluding %d ch (MCS %s)\n', ...
            sum(deadMask), mat2str(excluded));
    end
end

% -------------------------------------------------------------------------
function save_adj_heatmap(A, ttl, ~, layout, deadMask, outBase)
%SAVE_ADJ_HEATMAP Export adjacency CSV + JSON for Python heatmap plotting.
%   Excludes channels flagged in deadMask (reference + dead amplifier).
%   Produces Inkscape-ready PDFs via scripts/plot_heatmaps.py.

    [sortedLabels, sortOrder] = sort(layout.mcsLabels);

    % Map deadMask (unsorted) to sorted order
    keepSorted = ~deadMask(sortOrder);
    sortedLabels = sortedLabels(keepSorted);
    sortOrder    = sortOrder(keepSorted);
    A_sorted     = A(sortOrder, sortOrder);
    nCh          = numel(sortedLabels);

    % CSV: adjacency matrix with dead channels removed
    writematrix(A_sorted, [outBase '_data.csv']);

    % Tick labels: every electrode
    allLabels = arrayfun(@num2str, sortedLabels, 'UniformOutput', false);

    % JSON metadata
    symmetric = contains(ttl, '\Delta') || contains(ttl, 'Delta');
    meta = struct( ...
        'title',        ttl, ...
        'symmetric',    symmetric, ...
        'tick_pos',     1:nCh, ...
        'tick_labels',  {allLabels}, ...
        'n_channels',   nCh);
    fid = fopen([outBase '_meta.json'], 'w');
    fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
    fclose(fid);
end

% -------------------------------------------------------------------------
function save_network_panel(adjacency, nodeMetric, opt, ttl, edgeMode, ...
                            nodeCmap, showLabels, outBase)
    fig = create_panel_figure(6.0, 6.0);
    ax  = axes(fig);

    plot_network_on_mea(adjacency, ...
        'parent',           ax, ...
        'nodeMetric',       nodeMetric, ...
        'edgeMode',         edgeMode, ...
        'edgeThresholdPct', opt.edgeThresholdPct, ...
        'nodeCmap',         nodeCmap, ...
        'title',            ttl, ...
        'showLabels',       showLabels, ...
        'labelFontSize',    5);

    apply_nature_style(fig);
    save_figure(fig, outBase);
    close(fig);
end

