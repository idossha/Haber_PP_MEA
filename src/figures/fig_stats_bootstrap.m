function stats = fig_stats_bootstrap(study, varargin)
%FIG_STATS_BOOTSTRAP Hierarchical bootstrap on % change (median) per study.
%
%   stats = FIG_STATS_BOOTSTRAP(study) computes per-channel percent change
%   for both spike rate and burst rate using the same channel filters as
%   the corresponding violin scripts, runs a hierarchical bootstrap over
%   recordings then channels (default 10000 iterations), and returns a
%   struct with .spike and .burst fields:
%       .observed     - pooled median % change
%       .ci           - 95% bootstrap CI
%       .pApprox      - approximate two-sided p vs 0
%       .nDatasets    - number of recordings used
%       .nChannels    - total channel observations
%
%   Replaces the proof-of-concept figures_stats.m and prints a publication-
%   style summary line for each metric.
%
% INPUTS:
%   study  -  'doi' or 'ket'.
%
% Name-value options:
%   'pctMaxInclude'         -  default 1000
%   'excludeSilencedPct'    -  default false
%   'ignoreSilentChannels'  -  burst only, default true
%   'minRateThreshold'      -  burst only, default 0
%   'nBootstrap'            -  default 10000
%   'bootstrapSeed'         -  default 1
%   'channels'              -  default cfg.channels.default
%
% OUTPUTS:
%   stats  -  Struct with .spike and .burst fields (see above).

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pctMaxInclude',         cfg.pct.max_include);
    addParameter(p, 'excludeSilencedPct',    cfg.pct.exclude_silenced);
    addParameter(p, 'ignoreSilentChannels',  true);
    addParameter(p, 'minRateThreshold',      cfg.silent.min_rate_burst);
    addParameter(p, 'nBootstrap',            cfg.stats.n_bootstrap);
    addParameter(p, 'bootstrapSeed',         cfg.stats.bootstrap_seed);
    addParameter(p, 'channels',              cfg.channels.default);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    pairs = get_pairs_and_labels(cfg, study);

    nPairs = numel(pairs);
    spikePctByDataset = cell(nPairs, 1);
    burstPctByDataset = cell(nPairs, 1);

    for k = 1:nPairs
        [loadedB, loadedT] = load_pair_cache(pairs(k), cfg);

        chB = loadedB.channelsUsed(:)';
        chT = loadedT.channelsUsed(:)';
        [~, idxB] = ismember(opt.channels, chB);
        [~, idxT] = ismember(opt.channels, chT);
        valid = (idxB > 0) & (idxT > 0);
        if ~any(valid); continue; end

        rsB = loadedB.spikeRates(:);
        rsT = loadedT.spikeRates(:);

        % --- Spike pct ---
        bS = rsB(idxB(valid));
        tS = rsT(idxT(valid));
        ok = ~isnan(bS) & ~isnan(tS);
        bS = bS(ok); tS = tS(ok);
        keepS = ~(bS == 0 & tS == 0);
        bS = bS(keepS); tS = tS(keepS);
        if isempty(bS)
            pctS = [];
        else
            pctS = 100 * (tS - bS) ./ bS;
            infMask = isinf(pctS);
            if any(infMask)
                if ~isempty(opt.pctMaxInclude)
                    pctS(infMask) = opt.pctMaxInclude;
                else
                    pctS(infMask) = 1e6;
                end
            end
            if opt.excludeSilencedPct
                pctS = pctS(abs(pctS + 100) >= 1e-6);
            end
            if ~isempty(opt.pctMaxInclude)
                pctS = pctS(pctS <= opt.pctMaxInclude);
            end
            pctS = pctS(~isnan(pctS));
        end
        spikePctByDataset{k} = pctS;

        % --- Burst pct ---
        rbB = loadedB.burstRates(:);
        rbT = loadedT.burstRates(:);
        bB = rbB(idxB(valid));
        tB = rbT(idxT(valid));
        sB = rsB(idxB(valid));
        sT = rsT(idxT(valid));
        ok = ~isnan(bB) & ~isnan(tB) & ~isnan(sB) & ~isnan(sT);
        bB = bB(ok); tB = tB(ok); sB = sB(ok); sT = sT(ok);
        if opt.ignoreSilentChannels
            if opt.minRateThreshold == 0
                nonSilent = ~(bB == 0 & tB == 0);
            else
                nonSilent = ~(bB < opt.minRateThreshold & tB < opt.minRateThreshold);
            end
            bB = bB(nonSilent); tB = tB(nonSilent);
            sB = sB(nonSilent); sT = sT(nonSilent);
        end
        if isempty(bB)
            pctB = [];
        else
            keepSpike = ~((sB < 0 & sT < 0) | (sB == 0 & sT == 0));
            bB = bB(keepSpike); tB = tB(keepSpike);
            useB = bB > 0;
            bB = bB(useB); tB = tB(useB);
            if isempty(bB)
                pctB = [];
            else
                pctB = 100 * (tB - bB) ./ bB;
                if opt.excludeSilencedPct
                    pctB = pctB(abs(pctB + 100) >= 1e-6);
                end
                if ~isempty(opt.pctMaxInclude)
                    pctB = pctB(pctB <= opt.pctMaxInclude);
                end
                pctB = pctB(~isnan(pctB));
            end
        end
        burstPctByDataset{k} = pctB;
    end

    stats.spike = hierarchical_bootstrap_pct(spikePctByDataset, opt.nBootstrap, opt.bootstrapSeed);
    stats.burst = hierarchical_bootstrap_pct(burstPctByDataset, opt.nBootstrap, opt.bootstrapSeed);

    fprintf('\n--- fig_stats_bootstrap (study=%s, n=%d boots) ---\n', study, opt.nBootstrap);
    print_bootstrap_line('Spike rate', stats.spike);
    print_bootstrap_line('Burst rate', stats.burst);
    fprintf('--------------------------------------------------\n');

    % --- Numeric sidecar ------------------------------------------------
    out = struct( ...
        'study',                     study, ...
        'n_bootstrap',               opt.nBootstrap, ...
        'bootstrap_seed',            opt.bootstrapSeed, ...
        'spike_observed',            safe_field(stats.spike, 'observed'), ...
        'spike_ci_lo',               safe_ci(stats.spike, 1), ...
        'spike_ci_hi',               safe_ci(stats.spike, 2), ...
        'spike_p',                   safe_field(stats.spike, 'pApprox'), ...
        'spike_n_datasets',          safe_field(stats.spike, 'nDatasets'), ...
        'spike_n_channels',          safe_field(stats.spike, 'nChannels'), ...
        'burst_observed',            safe_field(stats.burst, 'observed'), ...
        'burst_ci_lo',               safe_ci(stats.burst, 1), ...
        'burst_ci_hi',               safe_ci(stats.burst, 2), ...
        'burst_p',                   safe_field(stats.burst, 'pApprox'), ...
        'burst_n_datasets',          safe_field(stats.burst, 'nDatasets'), ...
        'burst_n_channels',          safe_field(stats.burst, 'nChannels'));
    export_figure_stats(out, fullfile(cfg.paths.figures_out, ...
        sprintf('%s_stats_bootstrap', study)));
end

% -------------------------------------------------------------------------
function v = safe_field(s, f)
    if isfield(s, f); v = s.(f); else, v = NaN; end
end
function v = safe_ci(s, idx)
    if isfield(s, 'ci') && numel(s.ci) >= idx
        v = s.ci(idx);
    else
        v = NaN;
    end
end

% =========================================================================
function out = hierarchical_bootstrap_pct(pctByDataset, nBoot, rngSeed)
    pctByDataset = pctByDataset(cellfun(@(x) ~isempty(x), pctByDataset));
    for i = 1:numel(pctByDataset)
        pctByDataset{i} = pctByDataset{i}(:);
        pctByDataset{i} = pctByDataset{i}(~isnan(pctByDataset{i}));
    end
    pctByDataset = pctByDataset(~cellfun(@isempty, pctByDataset));

    if isempty(pctByDataset)
        out = struct('observed', NaN, 'ci', [NaN NaN], 'pApprox', NaN, ...
            'nDatasets', 0, 'nChannels', 0, 'bootStats', []);
        return;
    end

    rng(rngSeed);

    nD = numel(pctByDataset);
    nChannels = sum(cellfun(@numel, pctByDataset));
    obsAll = vertcat(pctByDataset{:});
    observed = median(obsAll, 'omitnan');

    bootStats = NaN(nBoot, 1);
    for b = 1:nBoot
        dsIdx = randi(nD, [nD, 1]);
        pooled = [];
        for j = 1:nD
            x = pctByDataset{dsIdx(j)};
            chIdx = randi(numel(x), [numel(x), 1]);
            pooled = [pooled; x(chIdx)]; %#ok<AGROW>
        end
        bootStats(b) = median(pooled, 'omitnan');
    end

    ci = prctile(bootStats, [2.5 97.5]);
    pApprox = 2 * min(mean(bootStats <= 0), mean(bootStats >= 0));
    pApprox = min(pApprox, 1);

    out = struct('observed', observed, 'ci', ci, 'pApprox', pApprox, ...
        'nDatasets', nD, 'nChannels', nChannels, 'bootStats', bootStats);
end

% -------------------------------------------------------------------------
function print_bootstrap_line(label, out)
    if out.nDatasets == 0
        fprintf('  %s: no valid datasets after filters.\n', label);
        return;
    end
    if out.pApprox < 0.001
        pStr = '<0.001';
    elseif out.pApprox < 0.01
        pStr = sprintf('%.3f', out.pApprox);
    else
        pStr = sprintf('%.3g', out.pApprox);
    end
    fprintf(['  %s: pooled median = %.4g%%, 95%% CI = [%.4g, %.4g]%%, ', ...
        'p approx = %s, n = %d datasets / %d channels\n'], ...
        label, out.observed, out.ci(1), out.ci(2), pStr, out.nDatasets, out.nChannels);
end
