%RUN_ENTROPY_ANALYSIS  Shannon entropy of per-channel binary spike trains.
%
%   For each channel, the spike train is binarized at 1 ms resolution (1 =
%   at least one spike in the bin, 0 = no spike). Shannon entropy is then
%   computed on the binary distribution:
%
%       H = -p0*log2(p0) - p1*log2(p1)
%
%   where p1 = fraction of bins containing a spike, p0 = 1 - p1.
%   Maximum entropy for a binary variable is 1 bit (when p0 = p1 = 0.5).
%   Silent channels (p1 = 0) have H = 0.
%
%   Additionally computes word entropy (Varley et al. 2024 style): the
%   spike train is divided into non-overlapping words of length L bins,
%   and the entropy of the word distribution is computed. This captures
%   temporal patterning beyond single-bin statistics.
%
%   Outputs:
%     figures_out/entropy_analysis.png / .pdf  — visualization
%     figures_out/entropy_stats.csv            — per-pair summary
%     Opens the figure for interactive assessment

cfg = project_config();
channels = cfg.channels.default;
binMs = 1;  % 1 ms bins (binary spike trains)
wordLength = 10;  % 10-ms words for word entropy (Varley uses similar)

[doiPairs, ~] = get_pairs_and_labels(cfg, 'doi');
[ketPairs, ~] = get_pairs_and_labels(cfg, 'ket');
nDoi = numel(doiPairs);
nKet = numel(ketPairs);

% =========================================================================
% Compute entropy for all conditions
% =========================================================================
fprintf('=== Shannon Entropy Analysis ===\n\n');

% Storage: per-channel entropy for each pair
doiBinH   = cell(nDoi, 2);  % {baseline, treatment} x nDoi
doiWordH  = cell(nDoi, 2);
ketBinH   = cell(nKet, 2);
ketWordH  = cell(nKet, 2);

fprintf('DOI arm:\n');
for k = 1:nDoi
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(doiPairs(k), channels, cfg);
    [doiBinH{k,1}, doiWordH{k,1}] = compute_channel_entropy(bSpikes, bMeta.durationSec, binMs, wordLength);
    [doiBinH{k,2}, doiWordH{k,2}] = compute_channel_entropy(tSpikes, tMeta.durationSec, binMs, wordLength);
    fprintf('  Pair %d/%d: baseline H=%.4f, DOI H=%.4f (median binary)\n', ...
        k, nDoi, median(doiBinH{k,1},'omitnan'), median(doiBinH{k,2},'omitnan'));
end

fprintf('\nKetanserin arm:\n');
for k = 1:nKet
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(ketPairs(k), channels, cfg);
    [ketBinH{k,1}, ketWordH{k,1}] = compute_channel_entropy(bSpikes, bMeta.durationSec, binMs, wordLength);
    [ketBinH{k,2}, ketWordH{k,2}] = compute_channel_entropy(tSpikes, tMeta.durationSec, binMs, wordLength);
    fprintf('  Pair %d/%d: baseline H=%.4f, Ket H=%.4f (median binary)\n', ...
        k, nKet, median(ketBinH{k,1},'omitnan'), median(ketBinH{k,2},'omitnan'));
end

% =========================================================================
% Compute per-well medians and paired stats
% =========================================================================
doiBinMedians  = zeros(nDoi, 2);
doiWordMedians = zeros(nDoi, 2);
ketBinMedians  = zeros(nKet, 2);
ketWordMedians = zeros(nKet, 2);

for k = 1:nDoi
    doiBinMedians(k,1)  = median(doiBinH{k,1}, 'omitnan');
    doiBinMedians(k,2)  = median(doiBinH{k,2}, 'omitnan');
    doiWordMedians(k,1) = median(doiWordH{k,1}, 'omitnan');
    doiWordMedians(k,2) = median(doiWordH{k,2}, 'omitnan');
end
for k = 1:nKet
    ketBinMedians(k,1)  = median(ketBinH{k,1}, 'omitnan');
    ketBinMedians(k,2)  = median(ketBinH{k,2}, 'omitnan');
    ketWordMedians(k,1) = median(ketWordH{k,1}, 'omitnan');
    ketWordMedians(k,2) = median(ketWordH{k,2}, 'omitnan');
end

% Print summary
fprintf('\n=== Summary ===\n');
fprintf('DOI Binary Entropy:  baseline %.4f -> DOI %.4f (delta %+.4f)\n', ...
    median(doiBinMedians(:,1)), median(doiBinMedians(:,2)), ...
    median(doiBinMedians(:,2) - doiBinMedians(:,1)));
fprintf('DOI Word Entropy:    baseline %.4f -> DOI %.4f (delta %+.4f)\n', ...
    median(doiWordMedians(:,1)), median(doiWordMedians(:,2)), ...
    median(doiWordMedians(:,2) - doiWordMedians(:,1)));
fprintf('Ket Binary Entropy:  baseline %.4f -> Ket %.4f (delta %+.4f)\n', ...
    median(ketBinMedians(:,1)), median(ketBinMedians(:,2)), ...
    median(ketBinMedians(:,2) - ketBinMedians(:,1)));
fprintf('Ket Word Entropy:    baseline %.4f -> Ket %.4f (delta %+.4f)\n', ...
    median(ketWordMedians(:,1)), median(ketWordMedians(:,2)), ...
    median(ketWordMedians(:,2) - ketWordMedians(:,1)));

% Wilcoxon signed-rank (per-well medians)
if nDoi >= 4
    [pDOIbin, ~, statsDOIbin]   = signrank(doiBinMedians(:,2), doiBinMedians(:,1));
    [pDOIword, ~, statsDOIword] = signrank(doiWordMedians(:,2), doiWordMedians(:,1));
else
    pDOIbin = NaN; pDOIword = NaN;
end
fprintf('\nDOI Wilcoxon p (binary): %.4f\n', pDOIbin);
fprintf('DOI Wilcoxon p (word):   %.4f\n', pDOIword);

% =========================================================================
% Write CSV
% =========================================================================
outDir = cfg.paths.figures_out;
rows = {};
for k = 1:nDoi
    rows{end+1,1} = 'doi'; rows{end,2} = k;
    rows{end,3} = doiBinMedians(k,1); rows{end,4} = doiBinMedians(k,2);
    rows{end,5} = doiWordMedians(k,1); rows{end,6} = doiWordMedians(k,2);
end
for k = 1:nKet
    rows{end+1,1} = 'ket'; rows{end,2} = k;
    rows{end,3} = ketBinMedians(k,1); rows{end,4} = ketBinMedians(k,2);
    rows{end,5} = ketWordMedians(k,1); rows{end,6} = ketWordMedians(k,2);
end
T = cell2table(rows, 'VariableNames', ...
    {'study','pair','bin_entropy_baseline','bin_entropy_treatment', ...
     'word_entropy_baseline','word_entropy_treatment'});
writetable(T, fullfile(outDir, 'entropy_stats.csv'));
fprintf('\nWrote entropy_stats.csv\n');

% =========================================================================
% Visualization: 2x2 panel figure
% =========================================================================
fig = figure('Position', [50 50 1100 900], 'Color', 'w');

cBase = [0.30 0.59 0.85];
cDOI  = [0.87 0.32 0.20];
cKet  = [0.30 0.69 0.47];

% --- Panel A: DOI binary entropy (paired dot plot) ---
subplot(2,2,1);
plot_paired(doiBinMedians, {'Baseline','DOI'}, cBase, cDOI);
ylabel('Binary entropy (bits)');
title('A. DOI — Binary Entropy', 'FontWeight', 'bold');
ylim([0 0.25]);

% --- Panel B: DOI word entropy (paired dot plot) ---
subplot(2,2,2);
plot_paired(doiWordMedians, {'Baseline','DOI'}, cBase, cDOI);
ylabel(sprintf('Word entropy (bits, L=%d ms)', wordLength));
title('B. DOI — Word Entropy', 'FontWeight', 'bold');

% --- Panel C: Ket binary entropy (paired dot plot) ---
subplot(2,2,3);
plot_paired(ketBinMedians, {'Baseline','DOI+Ket'}, cBase, cKet);
ylabel('Binary entropy (bits)');
title('C. DOI+Ket — Binary Entropy', 'FontWeight', 'bold');
ylim([0 0.25]);

% --- Panel D: Per-channel distributions (DOI, exemplar pair) ---
subplot(2,2,4);
exPair = 1;
bH = doiBinH{exPair, 1};
tH = doiBinH{exPair, 2};
% Remove silent channels
bH = bH(bH > 0);
tH = tH(tH > 0);
hold on;
histogram(bH, 20, 'FaceColor', cBase, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
histogram(tH, 20, 'FaceColor', cDOI,  'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;
xlabel('Binary entropy (bits)');
ylabel('Channel count');
legend('Baseline', 'DOI', 'Location', 'northeast');
title('D. Per-Channel Distribution (Pair 1)', 'FontWeight', 'bold');
set(gca, 'Box', 'off', 'TickDir', 'out');

% Save
exportgraphics(fig, fullfile(outDir, 'entropy_analysis.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'entropy_analysis.pdf'), 'ContentType', 'vector');
fprintf('Saved entropy_analysis.png and .pdf\n');
fprintf('\nDone. Figure is open for assessment.\n');

% =========================================================================
% Helper functions
% =========================================================================

function [binH, wordH] = compute_channel_entropy(spikesByChannel, durationSec, binMs, wordLen)
%COMPUTE_CHANNEL_ENTROPY  Binary and word entropy for each channel.
    nCh = numel(spikesByChannel);
    binSec = binMs / 1000;
    nBins = floor(durationSec / binSec);
    edges = (0:nBins) * binSec;

    binH  = zeros(nCh, 1);
    wordH = zeros(nCh, 1);

    for c = 1:nCh
        sp = spikesByChannel{c};
        if isempty(sp)
            binH(c)  = 0;
            wordH(c) = 0;
            continue;
        end

        % Binary spike train
        counts = histcounts(sp, edges);
        binary = double(counts > 0);

        % Binary entropy: H = -p0*log2(p0) - p1*log2(p1)
        p1 = mean(binary);
        p0 = 1 - p1;
        if p1 == 0 || p1 == 1
            binH(c) = 0;
        else
            binH(c) = -p0*log2(p0) - p1*log2(p1);
        end

        % Word entropy: divide into non-overlapping words of length wordLen
        nWords = floor(numel(binary) / wordLen);
        if nWords < 10
            wordH(c) = NaN;
            continue;
        end
        trimmed = binary(1:nWords*wordLen);
        words = reshape(trimmed, wordLen, nWords)';
        % Convert each word to a decimal index
        powers = 2.^(wordLen-1:-1:0);
        wordIdx = words * powers';
        % Count unique words
        [uniqueWords, ~, ic] = unique(wordIdx);
        nUnique = numel(uniqueWords);
        wordCounts = accumarray(ic, 1);
        probs = wordCounts / nWords;
        wordH(c) = -sum(probs .* log2(probs));
    end
end

function plot_paired(data, xlabels, c1, c2)
%PLOT_PAIRED  Paired dot plot with connecting lines.
    nPairs = size(data, 1);
    hold on;
    % Connecting lines
    for k = 1:nPairs
        if data(k,2) >= data(k,1)
            lc = [0.7 0.3 0.3]; % increase
        else
            lc = [0.3 0.3 0.7]; % decrease
        end
        plot([1 2], data(k,:), '-', 'Color', [lc 0.5], 'LineWidth', 1);
    end
    % Points
    scatter(ones(nPairs,1), data(:,1), 60, c1, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    scatter(2*ones(nPairs,1), data(:,2), 60, c2, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    % Medians
    plot(1, median(data(:,1)), 'd', 'MarkerSize', 12, 'MarkerFaceColor', c1, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    plot(2, median(data(:,2)), 'd', 'MarkerSize', 12, 'MarkerFaceColor', c2, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    hold off;
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', xlabels, 'FontSize', 11, ...
        'Box', 'off', 'TickDir', 'out');
end
