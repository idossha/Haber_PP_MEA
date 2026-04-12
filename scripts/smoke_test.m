function smoke_test()
%SMOKE_TEST Runtime sanity check for the Haber_PP_MEA active codebase.
%
%   Exercises the utilities, analysis modules, and stats helpers WITHOUT
%   touching any raw TDT data or the TDTMatlabSDK. Builds a synthetic v2.0
%   cache, writes it to a temp directory, and runs every loader /
%   connectivity / metric / stats function against it. Prints PASS/FAIL
%   lines and a final summary.
%
%   Usage:
%     matlab -nodisplay -nosplash -nodesktop -batch "addpath(genpath('src'));smoke_test"
%
% This test does NOT exercise: TDT loading, filtering, spike detection,
% or the figure scripts (those need a display and real spikes). It does
% exercise everything downstream of preprocess_and_save's output struct.

    tests = {};
    passCount = 0;
    failCount = 0;

    function record(name, ok, msg)
        if ok
            fprintf('  [PASS] %s\n', name);
            passCount = passCount + 1;
        else
            fprintf('  [FAIL] %s  --  %s\n', name, msg);
            failCount = failCount + 1;
        end
        tests{end+1} = struct('name', name, 'ok', ok, 'msg', msg); %#ok<AGROW>
    end

    fprintf('\n==== Haber_PP_MEA smoke test ====\n\n');

    % ---- Project config ------------------------------------------------
    try
        cfg = project_config();
        record('project_config() returns struct', isstruct(cfg), '');
        record('cfg.paths.root exists', isfield(cfg, 'paths') && isfield(cfg.paths, 'root'), '');
        record('cfg.datasets.doi has 6 pairs', ...
            numel(cfg.datasets.doi.baseline) == 6 && numel(cfg.datasets.doi.treatment) == 6, '');
        record('cfg.datasets.ket has 3 pairs', ...
            numel(cfg.datasets.ket.baseline) == 3 && numel(cfg.datasets.ket.treatment) == 3, '');
    catch ME
        record('project_config()', false, ME.message);
        cfg = struct();
    end

    % ---- MEA layout ----------------------------------------------------
    try
        layout = mea60_layout();
        record('mea60_layout() returns 60 channels', layout.nCh == 60, '');
        record('mea60_layout() grid is 8x8', isequal(layout.gridSize, [8 8]), '');
        record('distanceMatrix is 60x60', isequal(size(layout.distanceMatrix), [60 60]), '');
        record('distance matrix diagonal == 0', all(diag(layout.distanceMatrix) == 0), '');
    catch ME
        record('mea60_layout()', false, ME.message);
    end

    % ---- Pair datasets by #N suffix -----------------------------------
    try
        pairs = pair_datasets(cfg.datasets.doi.baseline, cfg.datasets.doi.treatment);
        record('pair_datasets DOI returns 6 pairs', numel(pairs) == 6, '');
        record('pairs(1).baseline ends in _#1', endsWith(pairs(1).baseline, '_#1'), '');
    catch ME
        record('pair_datasets', false, ME.message);
    end

    % ---- Build a synthetic v2.0 cache and exercise loaders -------------
    tmpCache = tempname;
    mkdir(tmpCache);
    cleanup = onCleanup(@() tryRmdir(tmpCache));

    try
        nCh = 60;
        fs  = 24414.0625;
        dur = 600;
        rng(42, 'twister');

        spikeTimes  = cell(1, nCh);
        spikeCounts = zeros(nCh, 1);
        spikeRates  = zeros(nCh, 1);
        burstStartIdx   = cell(1, nCh);
        burstEndIdx     = cell(1, nCh);
        burstStartTimes = cell(1, nCh);
        burstEndTimes   = cell(1, nCh);
        burstCounts = zeros(nCh, 1);
        burstRates  = zeros(nCh, 1);

        for c = 1:nCh
            % No Stats Toolbox -- use uniform sampling of spike times.
            lambda = 2 + 2 * rand();
            n = round(lambda * dur + sqrt(lambda * dur) * randn());
            n = max(n, 0);
            ts = sort(rand(n, 1) * dur);
            spikeTimes{c}  = ts;
            spikeCounts(c) = n;
            spikeRates(c)  = n * 60 / dur;

            % Dummy bursts: split the train into chunks of 5.
            if n >= 5
                nB = floor(n / 5);
                si = 1:5:(5*nB);
                ei = min(si + 4, n);
                burstStartIdx{c} = si;
                burstEndIdx{c}   = ei;
                burstStartTimes{c} = ts(si).';
                burstEndTimes{c}   = ts(ei).';
                burstCounts(c) = nB;
                burstRates(c)  = nB * 60 / dur;
            else
                burstStartIdx{c}   = zeros(1, 0);
                burstEndIdx{c}     = zeros(1, 0);
                burstStartTimes{c} = zeros(1, 0);
                burstEndTimes{c}   = zeros(1, 0);
            end
        end

        % Prototype cache struct.
        proto = struct();
        proto.version         = '2.0';
        proto.createdAt       = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
        proto.fs              = fs;
        proto.durationSec     = dur;
        proto.channelsUsed    = 1:nCh;
        proto.cfgUsed         = struct('synthetic', true);
        proto.spikeTimes      = spikeTimes;
        proto.spikeCounts     = spikeCounts;
        proto.spikeRates      = spikeRates;
        proto.burstStartIdx   = burstStartIdx;
        proto.burstEndIdx     = burstEndIdx;
        proto.burstStartTimes = burstStartTimes;
        proto.burstEndTimes   = burstEndTimes;
        proto.burstCounts     = burstCounts;
        proto.burstRates      = burstRates;

        % Fake cfg pointing at the tmp directory.
        fakeCfg = cfg;
        fakeCfg.paths.cache = tmpCache;

        % Write two caches for a one-pair "study" smoke test.
        b = 'IdoControl-230914-130200_#1';
        t = 'IdoDOI-230914-142502_#1';
        cB = proto; cB.datasetName = b; %#ok<NASGU>
        save(fullfile(tmpCache, cache_filename(b)), '-struct', 'cB', '-v7.3');
        cT = proto; cT.datasetName = t; %#ok<NASGU>
        save(fullfile(tmpCache, cache_filename(t)), '-struct', 'cT', '-v7.3');

        record('synthetic cache files written', ...
            exist(fullfile(tmpCache, cache_filename(b)), 'file') == 2 && ...
            exist(fullfile(tmpCache, cache_filename(t)), 'file') == 2, '');
    catch ME
        record('synthetic cache build', false, ME.message);
        return;
    end

    % ---- load_cache + load_pair_cache ---------------------------------
    try
        c1 = load_cache(b, fakeCfg);
        record('load_cache returns v2.0 struct', strcmp(c1.version, '2.0'), '');
        record('load_cache has spikeTimes cell', iscell(c1.spikeTimes) && numel(c1.spikeTimes) == nCh, '');
        record('load_cache spikeRates length', numel(c1.spikeRates) == nCh, '');

        pair.baseline  = b;
        pair.treatment = t;
        [bc, tc] = load_pair_cache(pair, fakeCfg);
        record('load_pair_cache returns two structs', ...
            isstruct(bc) && isstruct(tc), '');
    catch ME
        record('load_cache / load_pair_cache', false, ME.message);
    end

    % ---- load_pair_metric ---------------------------------------------
    try
        [bR, tR] = load_pair_metric(pair, 1:nCh, 'spikeRates', fakeCfg);
        record('load_pair_metric spikeRates returns vectors', ...
            numel(bR) == nCh && numel(tR) == nCh, '');
        [bR2, ~] = load_pair_metric(pair, 1:nCh, 'burstRates', fakeCfg);
        record('load_pair_metric burstRates works', numel(bR2) == nCh, '');
        try
            load_pair_metric(pair, 1:nCh, 'rates', fakeCfg); %#ok<NODEF>
            record('load_pair_metric rejects legacy "rates"', false, 'should have errored');
        catch
            record('load_pair_metric rejects legacy "rates"', true, '');
        end
    catch ME
        record('load_pair_metric', false, ME.message);
    end

    % ---- load_pair_spikes ---------------------------------------------
    try
        [bS, tS, bM, tM] = load_pair_spikes(pair, 1:nCh, fakeCfg);
        record('load_pair_spikes returns cells', ...
            iscell(bS) && iscell(tS) && numel(bS) == nCh, '');
        record('load_pair_spikes meta has datasetName', ...
            isfield(bM, 'datasetName') && isfield(tM, 'datasetName'), '');
    catch ME
        record('load_pair_spikes', false, ME.message);
    end

    % ---- transfer_entropy_d1 ------------------------------------------
    try
        teRes = transfer_entropy_d1(bS, dur, 'binMs', 3, 'minRateHz', 1);
        need = {'TE','meanTE','asymmetry','activeChannels','nActive', ...
                'binMs','minRateHz','topEdgesFrac'};
        missing = setdiff(need, fieldnames(teRes));
        record('transfer_entropy_d1 returns all required fields', ...
            isempty(missing), strjoin(missing, ','));
        record('transfer_entropy_d1 meanTE is finite non-negative scalar', ...
            isscalar(teRes.meanTE) && isfinite(teRes.meanTE) && teRes.meanTE >= 0, ...
            '');
    catch ME
        record('transfer_entropy_d1 basic', false, ME.message);
    end

    try
        % Two independent Poisson channels -> near-zero TE.
        rng(7, 'twister');
        indepCells = cell(1, 4);
        lam = 5;                                % 5 Hz
        for cc = 1:4
            nn = max(0, round(lam * dur + sqrt(lam * dur) * randn()));
            indepCells{cc} = sort(rand(nn, 1) * dur);
        end
        teIndep = transfer_entropy_d1(indepCells, dur, 'binMs', 3, 'minRateHz', 1);
        record('transfer_entropy_d1 near-zero on independent Poisson', ...
            isfinite(teIndep.meanTE) && teIndep.meanTE < 0.05, ...
            sprintf('meanTE=%.4f', teIndep.meanTE));
    catch ME
        record('transfer_entropy_d1 independence', false, ME.message);
    end

    % ---- burst_onset_coincidence --------------------------------------
    try
        boc = burst_onset_coincidence(burstStartTimes, 'windowMs', 50, 'minBursts', 3);
        need = {'C','windowMs','minBursts','channelsUsed','nPairs'};
        missing = setdiff(need, fieldnames(boc));
        record('burst_onset_coincidence returns required fields', ...
            isempty(missing), strjoin(missing, ','));
        record('burst_onset_coincidence C in [0,1]', ...
            isscalar(boc.C) && isfinite(boc.C) && boc.C >= 0 && boc.C <= 1, ...
            sprintf('C=%.4f', boc.C));
    catch ME
        record('burst_onset_coincidence basic', false, ME.message);
    end

    try
        % Identical burst trains across channels -> C ~= 1.0
        commonB = (10:20:580).';
        identCells = cell(1, 8);
        for cc = 1:8; identCells{cc} = commonB; end
        bocSame = burst_onset_coincidence(identCells, 'windowMs', 50, 'minBursts', 3);
        record('burst_onset_coincidence identical -> 1.0', ...
            isfinite(bocSame.C) && abs(bocSame.C - 1.0) < 0.05, ...
            sprintf('C=%.4f', bocSame.C));
    catch ME
        record('burst_onset_coincidence identical', false, ME.message);
    end

    try
        % Disjoint burst trains: each channel's bursts are shifted by 10 s
        % relative to its neighbour, so no two channels have any overlap
        % within the 50 ms window -> C should be ~0.
        disjCells = cell(1, 6);
        for cc = 1:6
            disjCells{cc} = (10 * cc : 80 : 480).';   % non-overlapping schedules
        end
        bocDisj = burst_onset_coincidence(disjCells, 'windowMs', 50, 'minBursts', 3);
        record('burst_onset_coincidence disjoint -> ~0', ...
            isfinite(bocDisj.C) && bocDisj.C < 0.1, ...
            sprintf('C=%.4f', bocDisj.C));
    catch ME
        record('burst_onset_coincidence disjoint', false, ME.message);
    end

    % ---- spike_entropy_shannon ----------------------------------------
    try
        entRes = spike_entropy_shannon(bS, dur);
        need = {'H','medianH','binMs','clipCounts'};
        missing = setdiff(need, fieldnames(entRes));
        record('spike_entropy_shannon returns required fields', ...
            isempty(missing) && numel(entRes.H) == nCh, ...
            strjoin(missing, ','));
        Hmax = log2(entRes.clipCounts + 1);
        record('spike_entropy_shannon medianH finite in [0, log2(alphabet)]', ...
            isfinite(entRes.medianH) && entRes.medianH >= 0 ...
                && entRes.medianH <= Hmax + 1e-9, ...
            sprintf('medianH=%.4f Hmax=%.4f', entRes.medianH, Hmax));
    catch ME
        record('spike_entropy_shannon', false, ME.message);
    end

    % ---- avalanches ---------------------------------------------------
    try
        avRes = avalanches(bS, dur);
        need = {'binMs','meanIeiMs','nAvalanches','sizeCounts','lifetimeCounts', ...
                'alphaLS','betaLS','alphaLS_rmse','betaLS_rmse','alphaMLE', ...
                'alphaMLE_xmin','alphaMLE_ksp','lrtPreferredModel','tailMassFrac', ...
                'sethnaResidual','classification'};
        missing = setdiff(need, fieldnames(avRes));
        record('avalanches returns all required fields', ...
            isempty(missing), strjoin(missing, ','));
        allowed = {'subcritical','critical','supercritical','indeterminate'};
        record('avalanches classification is one of allowed values', ...
            any(strcmp(avRes.classification, allowed)), ...
            sprintf('got %s', avRes.classification));
        record('avalanches alphaLS is finite negative scalar', ...
            isscalar(avRes.alphaLS) && isfinite(avRes.alphaLS) && avRes.alphaLS < 0, ...
            sprintf('alphaLS=%.4f', avRes.alphaLS));
        % Synthetic Poisson -> any valid classification; just check run.
        record('avalanches runs on synthetic Poisson without error', ...
            isstruct(avRes) && ~isempty(avRes.classification), '');
    catch ME
        record('avalanches', false, ME.message);
    end

    % ---- connectivity_xcorr -------------------------------------------
    try
        % Subsample to 12 channels for speed.
        sub = bS(1:12);
        conn = connectivity_xcorr(sub, 'durationSec', dur, 'binMs', 5, 'maxLagMs', 50);
        record('connectivity_xcorr returns 12x12 adjacency', ...
            isequal(size(conn.adjacency), [12 12]), '');
        record('connectivity_xcorr adjacency symmetric (values)', ...
            max(abs(conn.adjacency - conn.adjacency'), [], 'all', 'omitnan') < 1e-10, '');
        record('connectivity_xcorr diagonal is NaN', ...
            all(isnan(diag(conn.adjacency))), '');
    catch ME
        record('connectivity_xcorr', false, ME.message);
        conn.adjacency = rand(12); conn.adjacency = (conn.adjacency + conn.adjacency')/2;
    end

    % ---- network_metrics (extended) -----------------------------------
    try
        m = network_metrics(conn.adjacency, 'nullN', 5);
        need = {'density','clusteringMean','meanShortestPath', ...
                'globalEfficiency','modularity','smallWorldnessSigma', ...
                'nCommunities'};
        missing = setdiff(need, fieldnames(m));
        record('network_metrics returns all extended fields', isempty(missing), ...
            strjoin(missing, ','));
        record('network_metrics.density in [0,1]', ...
            m.density >= 0 && m.density <= 1, '');
    catch ME
        record('network_metrics', false, ME.message);
    end

    % ---- bh_fdr -------------------------------------------------------
    try
        p = [0.001 0.01 0.02 0.03 0.04 0.2 0.4];
        [adj, rej] = bh_fdr(p, 0.05);
        record('bh_fdr returns adjusted p values', ...
            numel(adj) == numel(p) && all(adj >= p - 1e-12), '');
        record('bh_fdr rejects some at q=0.05', any(rej), '');
    catch ME
        record('bh_fdr', false, ME.message);
    end

    % ---- hedges_g_av --------------------------------------------------
    try
        b1 = [1 2 3 4 5 6].';
        t1 = [2 3 4 5 6 7].';
        g = hedges_g_av(b1, t1);
        record('hedges_g_av returns finite scalar', isscalar(g) && isfinite(g), '');
    catch ME
        record('hedges_g_av', false, ME.message);
    end

    % ---- permutation_test_exact ---------------------------------------
    try
        a6 = [1.0 1.2 1.1 1.3 1.4 1.5];
        b3 = [0.2 0.3 0.25];
        res = permutation_test_exact(a6, b3);
        record('permutation_test_exact enumerates 84 perms', ...
            res.nPermutations == 84, sprintf('got %d', res.nPermutations));
        record('permutation_test_exact p in [0,1]', ...
            res.pValue >= 0 && res.pValue <= 1, '');
    catch ME
        record('permutation_test_exact', false, ME.message);
    end

    % ---- paired_stats -------------------------------------------------
    try
        b2 = (10:15).';
        t2 = b2 + (2:7).';
        ps = paired_stats(b2, t2, 'nBootstrap', 500);
        need = {'n','medianDelta','bootstrap','wilcoxon','hedgesGav'};
        missing = setdiff(need, fieldnames(ps));
        record('paired_stats returns all top-level fields', isempty(missing), ...
            strjoin(missing, ','));
        record('paired_stats bootstrap has ciDelta', ...
            isfield(ps.bootstrap, 'ciDelta') && numel(ps.bootstrap.ciDelta) == 2, '');
        record('paired_stats bootstrap CI is 2-element', ...
            ~isnan(ps.bootstrap.ciDelta(1)) || ~isnan(ps.bootstrap.ciDelta(2)), '');
        record('paired_stats hedgesGav finite', isfinite(ps.hedgesGav), '');
        if exist('signrank', 'file') == 2
            record('paired_stats wilcoxon.p in [0,1]', ...
                ~isnan(ps.wilcoxon.p) && ps.wilcoxon.p >= 0 && ps.wilcoxon.p <= 1, '');
        else
            record('paired_stats degrades gracefully w/o Stats Toolbox', ...
                isnan(ps.wilcoxon.p) && ~isempty(ps.wilcoxon.note), ps.wilcoxon.note);
        end
    catch ME
        record('paired_stats', false, ME.message);
    end

    % ---- Parse check: every active .m file must tokenize cleanly -----
    try
        srcDirs = { ...
            fullfile(cfg.paths.root, 'src', 'config'), ...
            fullfile(cfg.paths.root, 'src', 'utils'), ...
            fullfile(cfg.paths.root, 'src', 'pipeline'), ...
            fullfile(cfg.paths.root, 'src', 'figures'), ...
            fullfile(cfg.paths.root, 'src', 'analysis')};
        parseFailures = {};
        for d = 1:numel(srcDirs)
            files = dir(fullfile(srcDirs{d}, '*.m'));
            for f = 1:numel(files)
                fpath = fullfile(files(f).folder, files(f).name);
                try
                    mtree(fpath, '-file'); %#ok<NASGU>
                catch ME
                    parseFailures{end+1} = sprintf('%s: %s', ...
                        files(f).name, ME.message); %#ok<AGROW>
                end
            end
        end
        record('every .m file parses (mtree)', isempty(parseFailures), ...
            strjoin(parseFailures, '; '));
    catch ME
        record('parse check', false, ME.message);
    end

    % ---- Check toolbox availability -----------------------------------
    hasSignal = license('test', 'Signal_Toolbox') && exist('filtfilt', 'file') == 2;
    hasStats  = license('test', 'Statistics_Toolbox') && exist('signrank', 'file') == 2;
    record('Signal Processing Toolbox available', hasSignal, ...
        'required by preprocess_and_save and connectivity_xcorr');
    record('Statistics Toolbox available', hasStats, ...
        'required by paired_stats.signrank and fig_connectivity_summary.ksdensity');

    % ---- Exercise the end-to-end filter + detect path on synthetic signal.
    % This catches toolbox mismatches (e.g. designNotchPeakIIR vs iirnotch)
    % that only surface during a real preprocess_and_save run.
    if hasSignal
        try
            fsTest = 24414;
            duration = 2;
            tvec = (0:1/fsTest:duration - 1/fsTest)';
            x = 1e-4 * randn(numel(tvec), 1) + 0.05 * sin(2*pi*60*tvec);  % 60 Hz + noise

            % Notch stack construction (the line that broke the first run).
            stack = build_notch_filter(fsTest, cfg.filter.notch_f0, ...
                                               cfg.filter.notch_max, ...
                                               cfg.filter.notch_Q);
            record('build_notch_filter returns a stack', ...
                ~isempty(stack) && isfield(stack, 'B') && isfield(stack, 'A'), '');

            % Apply one notch with filtfilt.
            y = filtfilt(stack(1).B, stack(1).A, x);
            record('filtfilt notch runs on synthetic signal', ...
                numel(y) == numel(x) && all(isfinite(y)), '');

            % Bandpass construction and filtfilt.
            [bpB, bpA] = build_bandpass_filter(fsTest, ...
                cfg.filter.bp_low, cfg.filter.bp_high, cfg.filter.bp_order);
            yb = filtfilt(bpB, bpA, x);
            record('build_bandpass_filter + filtfilt runs', ...
                numel(yb) == numel(x) && all(isfinite(yb)), '');

            % Full notch-stack traversal (exact shape of load_and_filter).
            xFull = x;
            for sk = 1:numel(stack)
                xFull = filtfilt(stack(sk).B, stack(sk).A, xFull);
            end
            xFull = filtfilt(bpB, bpA, xFull);
            record('full filter chain (all notches + bandpass)', ...
                all(isfinite(xFull)), '');
        catch ME
            record('filter chain smoke', false, ...
                sprintf('%s: %s', ME.identifier, ME.message));
        end
    end

    % ---- Summary -------------------------------------------------------
    fprintf('\n==== Result: %d PASS, %d FAIL ====\n', passCount, failCount);
    if failCount > 0
        for k = 1:numel(tests)
            if ~tests{k}.ok
                fprintf('  FAIL: %s -- %s\n', tests{k}.name, tests{k}.msg);
            end
        end
        error('smoke_test:Failures', '%d smoke test(s) failed.', failCount);
    end
end

% =========================================================================
function tryRmdir(p)
    try
        if exist(p, 'dir'); rmdir(p, 's'); end
    catch
    end
end
