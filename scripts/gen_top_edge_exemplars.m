%GEN_TOP_EDGE_EXEMPLARS  Show top-N strongest baseline edges as correlograms.
%
%   For a given study ('doi' or 'ket') and pair index, compute the
%   cross-correlation adjacency matrix and plot the top N strongest
%   baseline edges as individual correlogram panels. Each panel is
%   annotated with the channel pair, correlation weight, inter-electrode
%   distance (pitch units), and peak lag (ms). This lets you manually
%   pick the best exemplar edge for the main-text figure.
%
%   Usage (from project root):
%       run scripts/gen_top_edge_exemplars.m
%   Or modify the parameters below and run.

%% ── Parameters (override these before calling) ─────────────────────────
if ~exist('study',   'var'); study   = 'ket'; end   % 'doi' or 'ket'
if ~exist('pairIdx', 'var'); pairIdx = 2;     end   % which pair
if ~exist('topN',    'var'); topN    = 10;    end   % number of top edges

%% ── Setup ───────────────────────────────────────────────────────────────
cfg = project_config();
[pairs, labels] = get_pairs_and_labels(cfg, study);
pair = pairs(pairIdx);

fprintf('Study: %s, Pair %d/%d\n', study, pairIdx, numel(pairs));
fprintf('  Baseline:  %s\n', pair.baseline);
fprintf('  Treatment: %s\n\n', pair.treatment);

channels = cfg.channels.default;
[bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, channels, cfg);

%% ── Cross-correlation ───────────────────────────────────────────────────
binMs    = cfg.connectivity.bin_ms;
maxLagMs = cfg.connectivity.max_lag_ms;
norm     = cfg.connectivity.normalization;

bRes = connectivity_xcorr(bSpikes, ...
    'durationSec', bMeta.durationSec, ...
    'binMs', binMs, 'maxLagMs', maxLagMs, 'normalization', norm);

tRes = connectivity_xcorr(tSpikes, ...
    'durationSec', tMeta.durationSec, ...
    'binMs', binMs, 'maxLagMs', maxLagMs, 'normalization', norm);

%% ── Electrode layout for distances ─────────────────────────────────────
layout = mea60_layout();
distMat = layout.distanceMatrix;

%% ── Find top N edges ────────────────────────────────────────────────────
A = bRes.adjacency;
nCh = size(A, 1);
mask = triu(true(nCh), 1);
vals = A;
vals(~mask) = -Inf;
vals(isnan(vals)) = -Inf;

[sortedVals, sortedIdx] = sort(vals(:), 'descend');
topN = min(topN, sum(isfinite(sortedVals)));

topEdges = zeros(topN, 2);
topWeights = zeros(topN, 1);
topDist = zeros(topN, 1);
topPeakLag = zeros(topN, 1);

for k = 1:topN
    [ii, jj] = ind2sub([nCh nCh], sortedIdx(k));
    topEdges(k, :) = [ii, jj];
    topWeights(k) = sortedVals(k);
    topDist(k) = distMat(ii, jj);
    topPeakLag(k) = bRes.peakLagMs(ii, jj);
end

%% ── Print summary table ────────────────────────────────────────────────
fprintf('%-5s  %-16s  %-10s  %-12s  %-10s\n', ...
    'Rank', 'Channels', 'Weight', 'Distance', 'PeakLag(ms)');
fprintf('%s\n', repmat('-', 1, 60));
for k = 1:topN
    ii = topEdges(k, 1); jj = topEdges(k, 2);
    flag = '';
    if topDist(k) <= 1.01
        flag = ' ** ADJACENT';
    end
    if abs(topPeakLag(k)) < 1.5
        flag = [flag ' * ZERO-LAG'];
    end
    fprintf('%-5d  ch %2d <-> ch %2d  %8.4f  %8.2f      %+6.1f%s\n', ...
        k, ii, jj, topWeights(k), topDist(k), topPeakLag(k), flag);
end

%% ── Plot correlograms ──────────────────────────────────────────────────
colors = paired_plot_colors(study);
nCols = min(5, topN);
nRows = ceil(topN / nCols);

fig = figure('Units', 'centimeters', 'Position', [2 2 nCols*8 nRows*6], ...
    'Color', 'w', 'Renderer', 'painters');
tl = tiledlayout(fig, nRows, nCols, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, sprintf('Top %d strongest baseline edges — %s pair %d', ...
    topN, upper(study), pairIdx), 'FontWeight', 'bold', 'FontSize', 12);

binSec = binMs / 1000;

for k = 1:topN
    ii = topEdges(k, 1);
    jj = topEdges(k, 2);

    ax = nexttile(tl);

    % Compute cross-correlograms
    [cB, lagsB] = local_pair_xcorr(bSpikes{ii}, bSpikes{jj}, binSec, ...
        bMeta.durationSec, maxLagMs);
    [cT, lagsT] = local_pair_xcorr(tSpikes{ii}, tSpikes{jj}, binSec, ...
        tMeta.durationSec, maxLagMs);

    hold(ax, 'on');
    plot(ax, lagsB, cB, 'Color', colors.baseline,  'LineWidth', 1.2);
    plot(ax, lagsT, cT, 'Color', colors.treatment, 'LineWidth', 1.2);
    xline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
    hold(ax, 'off');

    xlim(ax, [-maxLagMs maxLagMs]);
    xlabel(ax, 'Lag (ms)', 'FontSize', 7);
    ylabel(ax, 'Corr (z)', 'FontSize', 7);
    ax.FontSize = 7;

    % Title with info
    distStr = sprintf('d=%.1f', topDist(k));
    lagStr  = sprintf('lag=%+.0fms', topPeakLag(k));
    ttl = sprintf('#%d: ch%d\\leftrightarrowch%d  w=%.3f\n%s  %s', ...
        k, ii, jj, topWeights(k), distStr, lagStr);
    title(ax, ttl, 'FontSize', 7, 'FontWeight', 'normal');

    % Flag problematic edges
    if topDist(k) <= 1.01 || abs(topPeakLag(k)) < 1.5
        ax.Color = [1 0.95 0.95]; % light red background
    end

    if k == 1
        legend(ax, {labels.baseline, labels.treatment}, ...
            'Location', 'northeast', 'FontSize', 5, 'Box', 'off');
    end
end

%% ── Save ────────────────────────────────────────────────────────────────
outDir = output_path(cfg, study, 'connectivity', 'exemplars');
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, sprintf('%s_pair%d_top%d_edges.png', study, pairIdx, topN));
exportgraphics(fig, outFile, 'Resolution', 300);
fprintf('\nSaved: %s\n', outFile);

% Also save as PDF
outPdf = strrep(outFile, '.png', '.pdf');
exportgraphics(fig, outPdf, 'ContentType', 'vector');
fprintf('Saved: %s\n', outPdf);


%% ── Local helper ────────────────────────────────────────────────────────
function [c, lags] = local_pair_xcorr(tsI, tsJ, binSec, durationSec, maxLagMs)
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
