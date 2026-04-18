%PLOT_POPULATION_FIRING_RATE  Population firing rate over time with shaded SD.
%
%   Overlays three conditions on a single panel:
%     - Baseline (blue, pooled across DOI baseline wells)
%     - DOI (red/orange, pooled across DOI treatment wells)
%     - DOI + Ketanserin (green, pooled across Ket treatment wells)
%
%   For each well, all spike times across channels are pooled into a single
%   population spike train, then binned at 1 s resolution to produce a
%   population firing rate time series (Hz). The per-well traces are then
%   averaged across wells (mean +/- SD shaded ribbon).
%
%   Only the first 3 minutes of each recording are plotted.

cfg = project_config();
channels = cfg.channels.default;
plotDurationSec = 570;  % 9.5 min (longest common window across all 18 recordings)
binSec = 5;             % 5-second bins for smoother traces
smoothSpan = 5;         % additional moving-average smoothing (number of bins)

% =========================================================================
% Load spike times for all conditions
% =========================================================================
[doiPairs, ~]  = get_pairs_and_labels(cfg, 'doi');
[ketPairs, ~]  = get_pairs_and_labels(cfg, 'ket');

nDoi = numel(doiPairs);
nKet = numel(ketPairs);

nBins = floor(plotDurationSec / binSec);
edges = (0:nBins) * binSec;
tMid  = edges(1:end-1) + binSec/2;  % bin centers in seconds
tMin  = tMid / 60;                   % bin centers in minutes

% Collect population firing rates: each row = one well
baselineRates = zeros(nDoi, nBins);
doiRates      = zeros(nDoi, nBins);
ketRates      = zeros(nKet, nBins);

fprintf('Loading DOI pairs...\n');
for k = 1:nDoi
    [bSpikes, tSpikes, ~, ~] = load_pair_spikes(doiPairs(k), channels, cfg);
    baselineRates(k, :) = pool_and_bin(bSpikes, edges);
    doiRates(k, :)      = pool_and_bin(tSpikes, edges);
    fprintf('  Pair %d/%d done\n', k, nDoi);
end

fprintf('Loading Ketanserin pairs...\n');
for k = 1:nKet
    [~, tSpikes, ~, ~] = load_pair_spikes(ketPairs(k), channels, cfg);
    ketRates(k, :) = pool_and_bin(tSpikes, edges);
    fprintf('  Pair %d/%d done\n', k, nKet);
end

% =========================================================================
% Compute mean and SD across wells
% =========================================================================
bMean = movmean(mean(baselineRates, 1), smoothSpan);
bSEM  = movmean(std(baselineRates, 0, 1) / sqrt(nDoi), smoothSpan);

dMean = movmean(mean(doiRates, 1), smoothSpan);
dSEM  = movmean(std(doiRates, 0, 1) / sqrt(nDoi), smoothSpan);

kMean = movmean(mean(ketRates, 1), smoothSpan);
kSEM  = movmean(std(ketRates, 0, 1) / sqrt(nKet), smoothSpan);

% =========================================================================
% Plot
% =========================================================================
fig = figure('Position', [100 100 800 400], 'Color', 'w');
set(gca, 'Toolbar', []);
hold on;

% Colors
cBase = [0.30 0.59 0.85];   % blue
cDOI  = [0.87 0.32 0.20];   % red-orange
cKet  = [0.30 0.69 0.47];   % green

alpha = 0.25;

% Shaded ribbons (SD)
fill([tMin, fliplr(tMin)], ...
     [bMean + bSEM, fliplr(bMean - bSEM)], ...
     cBase, 'FaceAlpha', alpha, 'EdgeColor', 'none');
fill([tMin, fliplr(tMin)], ...
     [dMean + dSEM, fliplr(dMean - dSEM)], ...
     cDOI, 'FaceAlpha', alpha, 'EdgeColor', 'none');
fill([tMin, fliplr(tMin)], ...
     [kMean + kSEM, fliplr(kMean - kSEM)], ...
     cKet, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% Mean lines
h1 = plot(tMin, bMean, '-', 'Color', cBase, 'LineWidth', 1.5);
h2 = plot(tMin, dMean, '-', 'Color', cDOI,  'LineWidth', 1.5);
h3 = plot(tMin, kMean, '-', 'Color', cKet,  'LineWidth', 1.5);

hold off;

xlabel('Time (min)', 'FontSize', 12);
ylabel('Population firing rate (Hz)', 'FontSize', 12);
legend([h1, h2, h3], {'Baseline', 'DOI', 'DOI + Ket'}, ...
    'Location', 'northeast', 'FontSize', 10);
set(gca, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
xlim([0 plotDurationSec/60]);
ylim([0 inf]);

% Save
outDir = fullfile(cfg.paths.output, 'shared');
if ~exist(outDir, 'dir'); mkdir(outDir); end
exportgraphics(fig, fullfile(outDir, 'population_firing_rate.png'), ...
    'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'population_firing_rate.pdf'), ...
    'ContentType', 'vector');
fprintf('Saved to %s\n', outDir);

% =========================================================================
function rate = pool_and_bin(spikesByChannel, edges)
%POOL_AND_BIN Pool all spikes across channels and histogram into time bins.
    allSpikes = [];
    for c = 1:numel(spikesByChannel)
        sp = spikesByChannel{c};
        if ~isempty(sp)
            allSpikes = [allSpikes; sp(:)]; %#ok<AGROW>
        end
    end
    % Only count spikes within the plotting window
    maxT = edges(end);
    allSpikes = allSpikes(allSpikes <= maxT);
    counts = histcounts(allSpikes, edges);
    binWidth = edges(2) - edges(1);
    rate = counts / binWidth;  % spikes/sec = Hz
end
