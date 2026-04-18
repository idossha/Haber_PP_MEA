%RUN_SUPPLEMENTARY_CONNECTIVITY  Density sweep + STTC robustness analyses.
%
%   Runs two supplementary connectivity analyses for both DOI and ketanserin
%   arms, writes results to figures_out/ as CSV tables:
%
%   1. DENSITY SWEEP — recomputes all six graph-theoretic metrics at target
%      edge densities [5, 10, 15, 20, 25] % (Rubinov & Sporns 2010).
%
%   2. STTC — computes the spike-time tiling coefficient (Cutts & Eglen
%      2014) as a rate-independent connectivity estimator, then derives the
%      same six graph metrics and paired statistics.
%
%   Output files:
%     figures_out/<study>_density_sweep.csv
%     figures_out/<study>_sttc_metrics.csv
%     figures_out/<study>_sttc_vs_pearson.csv

cfg = project_config();
outDir = fullfile(cfg.paths.output, 'SI');
if ~exist(outDir, 'dir'); mkdir(outDir); end

studies = {'doi', 'ket'};

for s = 1:numel(studies)
    study = studies{s};
    fprintf('\n=== %s arm ===\n', upper(study));

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    nPairs = numel(pairs);
    channels = cfg.channels.default;

    % =====================================================================
    % PART 1: DENSITY SWEEP
    % =====================================================================
    fprintf('\n--- Density sweep ---\n');
    densities = cfg.connectivity.density_sweep;
    nDens = numel(densities);

    metricNames = {'density','clusteringMean','meanShortestPath', ...
                   'smallWorldnessSigma','globalEfficiency','modularity'};
    metricLabels = {'Edge density','Clustering coeff.','Char. path length', ...
                    'Small-worldness sigma','Global efficiency','Modularity Q'};

    % Collect: for each pair x condition x density, extract all 6 metrics.
    sweepRows = {};
    for k = 1:nPairs
        [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pairs(k), channels, cfg);
        fprintf('  Pair %d/%d: %s -> %s\n', k, nPairs, bMeta.datasetName, tMeta.datasetName);

        bRes = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, ...
            'binMs', cfg.connectivity.bin_ms, 'maxLagMs', cfg.connectivity.max_lag_ms);
        tRes = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, ...
            'binMs', cfg.connectivity.bin_ms, 'maxLagMs', cfg.connectivity.max_lag_ms);

        for d = 1:nDens
            targetDens = densities(d);
            bThresh = density_threshold(bRes.adjacency, targetDens);
            tThresh = density_threshold(tRes.adjacency, targetDens);

            bM = network_metrics(bRes.adjacency, 'threshold', bThresh, ...
                'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);
            tM = network_metrics(tRes.adjacency, 'threshold', tThresh, ...
                'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);

            for m = 1:numel(metricNames)
                mn = metricNames{m};
                sweepRows{end+1, 1} = study; %#ok<SAGROW>
                sweepRows{end, 2} = k;
                sweepRows{end, 3} = targetDens;
                sweepRows{end, 4} = metricLabels{m};
                sweepRows{end, 5} = bM.(mn);
                sweepRows{end, 6} = tM.(mn);
                sweepRows{end, 7} = tM.(mn) - bM.(mn);
            end
        end
    end

    sweepT = cell2table(sweepRows, 'VariableNames', ...
        {'study','pair','target_density','metric','baseline','treatment','delta'});
    sweepFile = fullfile(outDir, sprintf('%s_density_sweep.csv', study));
    writetable(sweepT, sweepFile);
    fprintf('  Wrote %s\n', sweepFile);

    % Print summary: for each density, show median delta for path length and modularity.
    fprintf('\n  Density sweep summary (median delta across pairs):\n');
    fprintf('  %8s  %18s  %14s\n', 'Density', 'Char. path length', 'Modularity Q');
    for d = 1:nDens
        mask_pl = sweepT.target_density == densities(d) & ...
                  strcmp(sweepT.metric, 'Char. path length');
        mask_mq = sweepT.target_density == densities(d) & ...
                  strcmp(sweepT.metric, 'Modularity Q');
        fprintf('  %8.0f%%  %+18.4f  %+14.4f\n', ...
            densities(d)*100, median(sweepT.delta(mask_pl)), median(sweepT.delta(mask_mq)));
    end

    % =====================================================================
    % PART 2: STTC CONNECTIVITY
    % =====================================================================
    fprintf('\n--- STTC analysis (delta = %d ms) ---\n', 50);

    sttcMetricRows = {};
    pearsonVals = struct();
    sttcVals    = struct();

    for mn_idx = 1:numel(metricNames)
        pearsonVals.(metricNames{mn_idx}).baseline  = zeros(nPairs, 1);
        pearsonVals.(metricNames{mn_idx}).treatment = zeros(nPairs, 1);
        sttcVals.(metricNames{mn_idx}).baseline     = zeros(nPairs, 1);
        sttcVals.(metricNames{mn_idx}).treatment    = zeros(nPairs, 1);
    end

    for k = 1:nPairs
        [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pairs(k), channels, cfg);
        fprintf('  Pair %d/%d: %s -> %s\n', k, nPairs, bMeta.datasetName, tMeta.datasetName);

        % --- Pearson (standard pipeline, for comparison) ---
        bPearson = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, ...
            'binMs', cfg.connectivity.bin_ms, 'maxLagMs', cfg.connectivity.max_lag_ms);
        tPearson = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, ...
            'binMs', cfg.connectivity.bin_ms, 'maxLagMs', cfg.connectivity.max_lag_ms);

        bPThresh = density_threshold(bPearson.adjacency, cfg.connectivity.edge_density);
        tPThresh = density_threshold(tPearson.adjacency, cfg.connectivity.edge_density);
        bPM = network_metrics(bPearson.adjacency, 'threshold', bPThresh, ...
            'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);
        tPM = network_metrics(tPearson.adjacency, 'threshold', tPThresh, ...
            'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);

        % --- STTC ---
        bSTTC = connectivity_sttc(bSpikes, 'durationSec', bMeta.durationSec, 'deltaMs', 50);
        tSTTC = connectivity_sttc(tSpikes, 'durationSec', tMeta.durationSec, 'deltaMs', 50);

        bSThresh = density_threshold(bSTTC.adjacency, cfg.connectivity.edge_density);
        tSThresh = density_threshold(tSTTC.adjacency, cfg.connectivity.edge_density);
        bSM = network_metrics(bSTTC.adjacency, 'threshold', bSThresh, ...
            'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);
        tSM = network_metrics(tSTTC.adjacency, 'threshold', tSThresh, ...
            'nullN', cfg.connectivity.null_N, 'nullSeed', cfg.connectivity.null_seed);

        for mn_idx = 1:numel(metricNames)
            mn = metricNames{mn_idx};
            pearsonVals.(mn).baseline(k)  = bPM.(mn);
            pearsonVals.(mn).treatment(k) = tPM.(mn);
            sttcVals.(mn).baseline(k)     = bSM.(mn);
            sttcVals.(mn).treatment(k)    = tSM.(mn);

            sttcMetricRows{end+1, 1} = study; %#ok<SAGROW>
            sttcMetricRows{end, 2} = k;
            sttcMetricRows{end, 3} = metricLabels{mn_idx};
            sttcMetricRows{end, 4} = bSM.(mn);
            sttcMetricRows{end, 5} = tSM.(mn);
            sttcMetricRows{end, 6} = tSM.(mn) - bSM.(mn);
        end
    end

    % Write STTC metrics table.
    sttcT = cell2table(sttcMetricRows, 'VariableNames', ...
        {'study','pair','metric','baseline','treatment','delta'});
    sttcFile = fullfile(outDir, sprintf('%s_sttc_metrics.csv', study));
    writetable(sttcT, sttcFile);
    fprintf('  Wrote %s\n', sttcFile);

    % Build Pearson vs STTC comparison table.
    compRows = {};
    fprintf('\n  STTC vs Pearson comparison (median delta, 10%% density):\n');
    fprintf('  %22s  %14s  %14s  %10s\n', 'Metric', 'Pearson delta', 'STTC delta', 'Same sign?');
    for mn_idx = 1:numel(metricNames)
        mn = metricNames{mn_idx};
        pDelta = pearsonVals.(mn).treatment - pearsonVals.(mn).baseline;
        sDelta = sttcVals.(mn).treatment    - sttcVals.(mn).baseline;
        medP = median(pDelta);
        medS = median(sDelta);
        sameSign = sign(medP) == sign(medS);
        fprintf('  %22s  %+14.4f  %+14.4f  %10s\n', ...
            metricLabels{mn_idx}, medP, medS, string(sameSign));

        for k = 1:nPairs
            compRows{end+1, 1} = study; %#ok<SAGROW>
            compRows{end, 2} = k;
            compRows{end, 3} = metricLabels{mn_idx};
            compRows{end, 4} = pDelta(k);
            compRows{end, 5} = sDelta(k);
        end
    end

    compT = cell2table(compRows, 'VariableNames', ...
        {'study','pair','metric','pearson_delta','sttc_delta'});
    compFile = fullfile(outDir, sprintf('%s_sttc_vs_pearson.csv', study));
    writetable(compT, compFile);
    fprintf('  Wrote %s\n', compFile);
end

fprintf('\n=== Done ===\n');

%% ========================================================================
function t = density_threshold(adjacency, targetDensity)
%DENSITY_THRESHOLD Compute threshold for target edge density.
    n = size(adjacency, 1);
    upperMask = triu(true(n), 1);
    vals = adjacency(upperMask);
    vals = vals(~isnan(vals));
    if isempty(vals)
        t = 0;
        return;
    end
    t = prctile(vals, 100 * (1 - targetDensity));
end
