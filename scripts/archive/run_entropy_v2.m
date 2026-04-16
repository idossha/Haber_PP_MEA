%RUN_ENTROPY_V2  Multi-scale entropy and complexity analysis.
%
%   Three approaches to capture temporal structure in spike trains:
%
%   1. Multi-scale Shannon entropy — binary spike trains at 10, 25, 50 ms
%      bins. Larger bins push p1 toward 0.5, making entropy more sensitive.
%
%   2. Lempel-Ziv complexity (LZC) — measures compressibility of the 1-ms
%      binary spike train. High LZC = random/complex; low LZC = repetitive/
%      patterned. Normalized to [0,1] relative to a random binary string of
%      equal length and density.
%
%   3. Sample entropy (SampEn) — measures self-similarity/predictability of
%      the instantaneous firing rate time series (spike counts in 50 ms
%      bins). Low SampEn = regular/predictable; high SampEn = irregular.
%      Standard parameters: m=2, r=0.2*SD (Richman & Moorman 2000).

cfg = project_config();
channels = cfg.channels.default;

[doiPairs, ~] = get_pairs_and_labels(cfg, 'doi');
[ketPairs, ~] = get_pairs_and_labels(cfg, 'ket');
nDoi = numel(doiPairs);
nKet = numel(ketPairs);

% Bin sizes for multi-scale Shannon entropy
binScales = [10 25 50];  % ms
nScales = numel(binScales);

% Storage
doiShannonMedians = zeros(nDoi, 2, nScales);  % (pair, condition, scale)
ketShannonMedians = zeros(nKet, 2, nScales);
doiLZCMedians     = zeros(nDoi, 2);
ketLZCMedians     = zeros(nKet, 2);
doiSampEnMedians  = zeros(nDoi, 2);
ketSampEnMedians  = zeros(nKet, 2);

fprintf('=== Multi-Scale Entropy & Complexity Analysis ===\n\n');

% =========================================================================
% DOI arm
% =========================================================================
fprintf('DOI arm:\n');
for k = 1:nDoi
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(doiPairs(k), channels, cfg);
    fprintf('  Pair %d/%d...', k, nDoi);

    for s = 1:nScales
        bH = channel_shannon(bSpikes, bMeta.durationSec, binScales(s));
        tH = channel_shannon(tSpikes, tMeta.durationSec, binScales(s));
        doiShannonMedians(k, 1, s) = median(bH, 'omitnan');
        doiShannonMedians(k, 2, s) = median(tH, 'omitnan');
    end

    bLZ = channel_lzc(bSpikes, bMeta.durationSec);
    tLZ = channel_lzc(tSpikes, tMeta.durationSec);
    doiLZCMedians(k, 1) = median(bLZ, 'omitnan');
    doiLZCMedians(k, 2) = median(tLZ, 'omitnan');

    bSE = channel_sampen(bSpikes, bMeta.durationSec);
    tSE = channel_sampen(tSpikes, tMeta.durationSec);
    doiSampEnMedians(k, 1) = median(bSE, 'omitnan');
    doiSampEnMedians(k, 2) = median(tSE, 'omitnan');

    fprintf(' done\n');
end

% =========================================================================
% KET arm
% =========================================================================
fprintf('\nKetanserin arm:\n');
for k = 1:nKet
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(ketPairs(k), channels, cfg);
    fprintf('  Pair %d/%d...', k, nKet);

    for s = 1:nScales
        bH = channel_shannon(bSpikes, bMeta.durationSec, binScales(s));
        tH = channel_shannon(tSpikes, tMeta.durationSec, binScales(s));
        ketShannonMedians(k, 1, s) = median(bH, 'omitnan');
        ketShannonMedians(k, 2, s) = median(tH, 'omitnan');
    end

    bLZ = channel_lzc(bSpikes, bMeta.durationSec);
    tLZ = channel_lzc(tSpikes, tMeta.durationSec);
    ketLZCMedians(k, 1) = median(bLZ, 'omitnan');
    ketLZCMedians(k, 2) = median(tLZ, 'omitnan');

    bSE = channel_sampen(bSpikes, bMeta.durationSec);
    tSE = channel_sampen(tSpikes, tMeta.durationSec);
    ketSampEnMedians(k, 1) = median(bSE, 'omitnan');
    ketSampEnMedians(k, 2) = median(tSE, 'omitnan');

    fprintf(' done\n');
end

% =========================================================================
% Summary statistics
% =========================================================================
fprintf('\n=== Summary ===\n\n');

fprintf('--- Multi-scale Shannon entropy (median across wells) ---\n');
for s = 1:nScales
    bm = doiShannonMedians(:,1,s); tm = doiShannonMedians(:,2,s);
    d = tm - bm;
    if nDoi >= 4
        pval = signrank(tm, bm);
    else
        pval = NaN;
    end
    fprintf('  DOI %3d ms bin:  baseline %.4f -> DOI %.4f  delta %+.4f  p=%.4f  (%d/%d up)\n', ...
        binScales(s), median(bm), median(tm), median(d), pval, sum(d>0), nDoi);
end
for s = 1:nScales
    bm = ketShannonMedians(:,1,s); tm = ketShannonMedians(:,2,s);
    d = tm - bm;
    fprintf('  Ket %3d ms bin:  baseline %.4f -> Ket %.4f  delta %+.4f  (%d/%d up)\n', ...
        binScales(s), median(bm), median(tm), median(d), sum(d>0), nKet);
end

fprintf('\n--- Lempel-Ziv Complexity (normalized) ---\n');
bm = doiLZCMedians(:,1); tm = doiLZCMedians(:,2); d = tm - bm;
pLZC = signrank(tm, bm);
fprintf('  DOI:  baseline %.4f -> DOI %.4f  delta %+.4f  p=%.4f  (%d/%d up)\n', ...
    median(bm), median(tm), median(d), pLZC, sum(d>0), nDoi);
bm = ketLZCMedians(:,1); tm = ketLZCMedians(:,2); d = tm - bm;
fprintf('  Ket:  baseline %.4f -> Ket %.4f  delta %+.4f  (%d/%d up)\n', ...
    median(bm), median(tm), median(d), sum(d>0), nKet);

fprintf('\n--- Sample Entropy (m=2, r=0.2*SD, 50ms bins) ---\n');
bm = doiSampEnMedians(:,1); tm = doiSampEnMedians(:,2); d = tm - bm;
pSE = signrank(tm, bm);
fprintf('  DOI:  baseline %.4f -> DOI %.4f  delta %+.4f  p=%.4f  (%d/%d up)\n', ...
    median(bm), median(tm), median(d), pSE, sum(d>0), nDoi);
bm = ketSampEnMedians(:,1); tm = ketSampEnMedians(:,2); d = tm - bm;
fprintf('  Ket:  baseline %.4f -> Ket %.4f  delta %+.4f  (%d/%d up)\n', ...
    median(bm), median(tm), median(d), sum(d>0), nKet);

% =========================================================================
% Visualization: 3x2 panel (DOI left, KET right)
% =========================================================================
fig = figure('Position', [50 50 1200 900], 'Color', 'w');

cBase = [0.30 0.59 0.85];
cDOI  = [0.87 0.32 0.20];
cKet  = [0.30 0.69 0.47];

% --- Row 1: Multi-scale Shannon entropy (best scale = 50 ms) ---
subplot(3,2,1);
plot_paired(doiShannonMedians(:,:,3), {'Baseline','DOI'}, cBase, cDOI);
ylabel('Shannon H (bits, 50 ms)');
title('A. DOI — Shannon Entropy (50 ms)', 'FontWeight', 'bold');

subplot(3,2,2);
plot_paired(ketShannonMedians(:,:,3), {'Baseline','DOI+Ket'}, cBase, cKet);
ylabel('Shannon H (bits, 50 ms)');
title('B. DOI+Ket — Shannon Entropy (50 ms)', 'FontWeight', 'bold');

% --- Row 2: Lempel-Ziv Complexity ---
subplot(3,2,3);
plot_paired(doiLZCMedians, {'Baseline','DOI'}, cBase, cDOI);
ylabel('Normalized LZC');
title('C. DOI — Lempel-Ziv Complexity', 'FontWeight', 'bold');

subplot(3,2,4);
plot_paired(ketLZCMedians, {'Baseline','DOI+Ket'}, cBase, cKet);
ylabel('Normalized LZC');
title('D. DOI+Ket — Lempel-Ziv Complexity', 'FontWeight', 'bold');

% --- Row 3: Sample Entropy ---
subplot(3,2,5);
plot_paired(doiSampEnMedians, {'Baseline','DOI'}, cBase, cDOI);
ylabel('SampEn');
title('E. DOI — Sample Entropy', 'FontWeight', 'bold');

subplot(3,2,6);
plot_paired(ketSampEnMedians, {'Baseline','DOI+Ket'}, cBase, cKet);
ylabel('SampEn');
title('F. DOI+Ket — Sample Entropy', 'FontWeight', 'bold');

% Save
outDir = cfg.paths.figures_out;
exportgraphics(fig, fullfile(outDir, 'entropy_v2_analysis.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'entropy_v2_analysis.pdf'), 'ContentType', 'vector');
fprintf('\nSaved entropy_v2_analysis.png and .pdf\n');

% =========================================================================
% Helper: Shannon entropy at a given bin size
% =========================================================================
function H = channel_shannon(spikesByChannel, durationSec, binMs)
    nCh = numel(spikesByChannel);
    binSec = binMs / 1000;
    nBins = floor(durationSec / binSec);
    edges = (0:nBins) * binSec;
    H = zeros(nCh, 1);
    for c = 1:nCh
        sp = spikesByChannel{c};
        if isempty(sp); H(c) = 0; continue; end
        counts = histcounts(sp, edges);
        binary = double(counts > 0);
        p1 = mean(binary);
        p0 = 1 - p1;
        if p1 == 0 || p1 == 1
            H(c) = 0;
        else
            H(c) = -p0*log2(p0) - p1*log2(p1);
        end
    end
end

% =========================================================================
% Helper: Lempel-Ziv complexity (normalized)
% =========================================================================
function C = channel_lzc(spikesByChannel, durationSec)
    nCh = numel(spikesByChannel);
    binSec = 0.001; % 1 ms
    nBins = floor(durationSec / binSec);
    edges = (0:nBins) * binSec;
    C = zeros(nCh, 1);
    for c = 1:nCh
        sp = spikesByChannel{c};
        if isempty(sp) || numel(sp) < 2
            C(c) = 0;
            continue;
        end
        counts = histcounts(sp, edges);
        binary = counts > 0;
        C(c) = lempel_ziv(binary);
    end
end

% =========================================================================
% Lempel-Ziv complexity (Kaspar & Schuster 1987, normalized)
% =========================================================================
function c = lempel_ziv(s)
    n = numel(s);
    if n == 0; c = 0; return; end
    % Count the number of distinct words
    nWords = 1;
    i = 1;
    k = 1;
    kmax = 1;
    while i + k <= n
        if s(i + k) ~= s(nWords + k)  % mismatch in extension
            if k > kmax
                kmax = k;
            end
            i = i + 1;
            if i == nWords + 1
                nWords = nWords + kmax + 1;
                i = nWords;
                k = 1;
                kmax = 1;
            else
                k = 1;
            end
        else
            k = k + 1;
        end
    end
    nWords = nWords + 1;
    % Normalize by theoretical upper bound for random binary string
    if n > 1
        b = n / log2(n);
        c = nWords / b;
    else
        c = 1;
    end
end

% =========================================================================
% Helper: Sample entropy (Richman & Moorman 2000)
% =========================================================================
function SE = channel_sampen(spikesByChannel, durationSec)
    nCh = numel(spikesByChannel);
    binMs = 200; % 200 ms bins — keeps N manageable for O(N^2)
    binSec = binMs / 1000;
    nBins = min(floor(durationSec / binSec), 1500); % cap at 1500 bins (5 min)
    edges = (0:nBins) * binSec;
    m = 2;      % embedding dimension
    SE = zeros(nCh, 1);
    for c = 1:nCh
        sp = spikesByChannel{c};
        if isempty(sp) || numel(sp) < 10
            SE(c) = NaN;
            continue;
        end
        counts = histcounts(sp, edges);
        rateSeries = double(counts);
        SE(c) = sampen(rateSeries, m);
    end
end

% =========================================================================
% Sample entropy implementation
% =========================================================================
function e = sampen(x, m)
    N = numel(x);
    if N < m + 2
        e = NaN;
        return;
    end
    r = 0.2 * std(x);
    if r == 0
        e = NaN;
        return;
    end
    % Count template matches for dimension m and m+1
    Bcount = 0; % m-length matches
    Acount = 0; % (m+1)-length matches
    for i = 1:N-m
        for j = i+1:N-m
            % Check m-length match
            if max(abs(x(i:i+m-1) - x(j:j+m-1))) < r
                Bcount = Bcount + 1;
                % Check (m+1)-length match
                if i+m <= N && j+m <= N
                    if abs(x(i+m) - x(j+m)) < r
                        Acount = Acount + 1;
                    end
                end
            end
        end
    end
    if Bcount == 0
        e = NaN;
    else
        e = -log(Acount / Bcount);
    end
end

% =========================================================================
% Paired dot plot helper
% =========================================================================
function plot_paired(data, xlabels, c1, c2)
    nPairs = size(data, 1);
    hold on;
    for k = 1:nPairs
        if data(k,2) >= data(k,1)
            lc = [0.7 0.3 0.3];
        else
            lc = [0.3 0.3 0.7];
        end
        plot([1 2], data(k,:), '-', 'Color', [lc 0.5], 'LineWidth', 1);
    end
    scatter(ones(nPairs,1), data(:,1), 60, c1, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    scatter(2*ones(nPairs,1), data(:,2), 60, c2, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    plot(1, median(data(:,1),'omitnan'), 'd', 'MarkerSize', 12, 'MarkerFaceColor', c1, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    plot(2, median(data(:,2),'omitnan'), 'd', 'MarkerSize', 12, 'MarkerFaceColor', c2, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hold off;
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', xlabels, 'FontSize', 11, ...
        'Box', 'off', 'TickDir', 'out');
end
