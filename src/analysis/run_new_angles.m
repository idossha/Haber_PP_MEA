function results = run_new_angles(study, varargin)
%RUN_NEW_ANGLES Per-study driver for the four Figure 6 "new angle" metrics.
%
%   results = RUN_NEW_ANGLES(study) computes, for every baseline/treatment
%   pair in STUDY, the four new-angle metrics we promote into Figure 6:
%       1. Delay-1 transfer entropy (Ito 2011 p.4 D1TE)
%       2. Burst-onset cross-channel coincidence (Brofiga 2023;
%          Chiappalone 2006 burst-propagation convention)
%       3. Per-channel Shannon entropy of binned spike counts
%          (Varley 2024 §2.2 50 ms bin, clip-3)
%       4. Neuronal-avalanche size/lifetime statistics (Pasquale 2008
%          p.1357 LS + Massobrio 2015 p.14 Clauset MLE)
%
%   Results are cached to <cfg.paths.figures_out>/<study>_new_angles.mat
%   so the five Figure 6 panel scripts (fig_new_angles_A..E) are I/O-bound
%   instead of recomputing TE (~10 s / recording) on every render.
%
%   results = RUN_NEW_ANGLES('doi')      runs the DOI arm (n = 6 pairs)
%   results = RUN_NEW_ANGLES('ket')      runs the ketanserin arm (n = 3 pairs)
%
% Name-value options:
%   'forceRecompute'  - recompute even if cache exists (default false)
%   'verbose'         - print per-pair progress (default true)
%
% INPUTS:
%   study  -  'doi' or 'ket' (case-insensitive).
%
% OUTPUTS:
%   results.study                 'doi' | 'ket'
%   results.pairLabels            1 x nPairs cell of "pair N: <basename>"
%   results.nPairs                scalar
%   results.te.baseline           nPairs x 1 mean TE at baseline (bits)
%   results.te.treatment          nPairs x 1 mean TE at treatment (bits)
%   results.te.asymBaseline       nPairs x 1 asymmetry at baseline
%   results.te.asymTreatment      nPairs x 1 asymmetry at treatment
%   results.burstSync.baseline    nPairs x 1 C at baseline
%   results.burstSync.treatment   nPairs x 1 C at treatment
%   results.entropy.baseline      nPairs x 1 median H at baseline (bits)
%   results.entropy.treatment     nPairs x 1 median H at treatment (bits)
%   results.avalanches.baseline   nPairs x 1 struct with fields alphaLS,
%                                  betaLS, alphaMLE, alphaMLE_ksp,
%                                  lrtPreferredModel, classification,
%                                  sizeCounts, lifetimeCounts,
%                                  sethnaResidual, binMs, meanIeiMs,
%                                  nAvalanches, tailMassFrac
%   results.avalanches.treatment  nPairs x 1 struct (same shape)
%   results.durationSec.baseline  nPairs x 1 (seconds)
%   results.durationSec.treatment nPairs x 1 (seconds)
%
% References:
%   Ito S, Hansen ME, Heiland R, Lumsdaine A, Litke AM, Beggs JM.
%     Extending transfer entropy improves identification of effective
%     connectivity in a spiking cortical network model. PLoS ONE
%     2011;6:e27431. (§p.4 D1TE definition.)
%   Varley TF et al. Information processing dynamics in neural networks
%     of macroscopic brain organoids. J Neural Eng 2024. (§2.2 Shannon H;
%     §2.3 multivariate TE.)
%   Pasquale V, Massobrio P, Bologna LL, Chiappalone M, Martinoia S.
%     Self-organization and neuronal avalanches in networks of dissociated
%     cortical neurons. Neuroscience 2008;153:1354.
%   Massobrio P, Pasquale V, Martinoia S. Self-organized criticality in
%     cortical assemblies. Sci Rep 2015;5:10578.
%   Brofiga M et al. On the functional role of excitatory and inhibitory
%     populations on the overall network dynamics. 2023.
%   Chiappalone M et al. Dissociated cortical networks show spontaneously
%     correlated activity patterns during in vitro development. Brain
%     Res 2006;1093:41-53.
%
% See also: TRANSFER_ENTROPY_D1, BURST_ONSET_COINCIDENCE,
%           SPIKE_ENTROPY_SHANNON, AVALANCHES, PAIRED_STATS.

    %%% ----- parse inputs --------------------------------------------------
    cfg = project_config();

    ip = inputParser;
    addRequired(ip,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(ip, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'verbose',        true,  @(x) islogical(x) && isscalar(x));
    parse(ip, study, varargin{:});
    opt = ip.Results;
    study = lower(opt.study);

    %%% ----- cache path ---------------------------------------------------
    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end
    cachePath = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles.mat', study));

    %%% ----- load from cache if available ---------------------------------
    if ~opt.forceRecompute && exist(cachePath, 'file') == 2
        try
            cached = load(cachePath);
            if isfield(cached, 'study') && strcmpi(cached.study, study)
                if opt.verbose
                    fprintf('run_new_angles(%s): loaded cache %s\n', ...
                        study, cachePath);
                end
                results = cached;
                return;
            end
        catch ME
            warning('run_new_angles:CacheLoadFailed', ...
                'Failed to load cache %s (%s). Recomputing.', ...
                cachePath, ME.message);
        end
    end

    %%% ----- compute per-pair metrics -------------------------------------
    [pairs, ~] = get_pairs_and_labels(cfg, study);
    nPairs = numel(pairs);
    channels = cfg.channels.default;

    pairLabels       = cell(1, nPairs);
    teBase           = nan(nPairs, 1);
    teTreat          = nan(nPairs, 1);
    asymBase         = nan(nPairs, 1);
    asymTreat        = nan(nPairs, 1);
    bsBase           = nan(nPairs, 1);
    bsTreat          = nan(nPairs, 1);
    entBase          = nan(nPairs, 1);
    entTreat         = nan(nPairs, 1);
    durBase          = nan(nPairs, 1);
    durTreat         = nan(nPairs, 1);
    protoAv          = compact_av_row(empty_av_struct());
    avBase           = repmat(protoAv, nPairs, 1);
    avTreat          = repmat(protoAv, nPairs, 1);

    for pp = 1:nPairs
        tStart = tic;
        pairLabels{pp} = sprintf('pair %d: %s', pp, pairs(pp).baseline);

        % --- load spikes + caches (baseline and treatment) --------------
        [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pairs(pp), channels, cfg);
        [bCache, tCache] = load_pair_cache(pairs(pp), cfg);

        durBase(pp)  = bMeta.durationSec;
        durTreat(pp) = tMeta.durationSec;

        % --- TE (Ito 2011 D1TE) ----------------------------------------
        teB = transfer_entropy_d1(bSpikes, bMeta.durationSec);
        teT = transfer_entropy_d1(tSpikes, tMeta.durationSec);
        teBase(pp)    = teB.meanTE;
        teTreat(pp)   = teT.meanTE;
        asymBase(pp)  = teB.asymmetry;
        asymTreat(pp) = teT.asymmetry;

        % --- Burst-onset coincidence (Brofiga 2023 / Chiappalone 2006) --
        bstB = reindex_bursts(bCache, channels);
        bstT = reindex_bursts(tCache, channels);
        bcB = burst_onset_coincidence(bstB);
        bcT = burst_onset_coincidence(bstT);
        bsBase(pp)  = bcB.C;
        bsTreat(pp) = bcT.C;

        % --- Shannon entropy (Varley 2024 §2.2) -------------------------
        eB = spike_entropy_shannon(bSpikes, bMeta.durationSec);
        eT = spike_entropy_shannon(tSpikes, tMeta.durationSec);
        entBase(pp)  = eB.medianH;
        entTreat(pp) = eT.medianH;

        % --- Avalanches (primary bin only; SI runs sweep separately) ----
        avB = avalanches(bSpikes, bMeta.durationSec, 'runSweep', false);
        avT = avalanches(tSpikes, tMeta.durationSec, 'runSweep', false);
        avBase(pp,1)  = compact_av_row(avB);
        avTreat(pp,1) = compact_av_row(avT);

        if opt.verbose
            fprintf(['[pair %d/%d] base TE=%.4f C=%.3f H=%.3f aLS=%+.2f | ', ...
                     'treat TE=%.4f C=%.3f H=%.3f aLS=%+.2f  (t=%.1fs)\n'], ...
                pp, nPairs, ...
                teBase(pp), bsBase(pp), entBase(pp), avBase(pp).alphaLS, ...
                teTreat(pp), bsTreat(pp), entTreat(pp), avTreat(pp).alphaLS, ...
                toc(tStart));
        end
    end

    %%% ----- pack results ------------------------------------------------
    results = struct();
    results.version     = '1.0';
    results.createdAt   = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    results.study       = study;
    results.pairLabels  = pairLabels;
    results.nPairs      = nPairs;
    results.te          = struct( ...
        'baseline',     teBase, ...
        'treatment',    teTreat, ...
        'asymBaseline', asymBase, ...
        'asymTreatment',asymTreat);
    results.burstSync   = struct( ...
        'baseline',     bsBase, ...
        'treatment',    bsTreat);
    results.entropy     = struct( ...
        'baseline',     entBase, ...
        'treatment',    entTreat);
    results.avalanches  = struct( ...
        'baseline',     avBase, ...
        'treatment',    avTreat);
    results.durationSec = struct( ...
        'baseline',     durBase, ...
        'treatment',    durTreat);

    %%% ----- save cache --------------------------------------------------
    try
        save(cachePath, '-struct', 'results', '-v7.3');
        if opt.verbose
            fprintf('run_new_angles(%s): saved cache %s\n', study, cachePath);
        end
    catch ME
        warning('run_new_angles:CacheSaveFailed', ...
            'Failed to save cache %s (%s).', cachePath, ME.message);
    end
end

% =========================================================================
function bst = reindex_bursts(cacheStruct, channels)
%REINDEX_BURSTS Align burstStartTimes to the requested channel list.
%
%   The cache stores burstStartTimes as a 1xnCh_used cell aligned to
%   cacheStruct.channelsUsed. We must map that to the caller's
%   `channels` ordering so burst_onset_coincidence sees one slot per
%   requested channel (empty when missing), matching the TE and entropy
%   calls above which already pass through load_pair_spikes.
    nCh = numel(channels);
    bst = repmat({zeros(0, 1)}, 1, nCh);
    chUsed = cacheStruct.channelsUsed(:)';
    [~, idxC] = ismember(channels, chUsed);
    src = cacheStruct.burstStartTimes;
    for k = 1:nCh
        if idxC(k) > 0 && idxC(k) <= numel(src)
            v = src{idxC(k)};
            bst{k} = v(:);
        end
    end
end

% =========================================================================
function r = compact_av_row(av)
%COMPACT_AV_ROW Distill avalanches() output to the fields we cache.
%
%   We keep only the scalars + small histograms that fig_new_angles_A
%   and fig_new_angles_E consume, so the cache file stays compact.
    r.alphaLS           = get_field(av, 'alphaLS',           NaN);
    r.betaLS            = get_field(av, 'betaLS',            NaN);
    r.alphaMLE          = get_field(av, 'alphaMLE',          NaN);
    r.alphaMLE_ksp      = get_field(av, 'alphaMLE_ksp',      NaN);
    r.lrtPreferredModel = get_field(av, 'lrtPreferredModel', 'power_law');
    r.classification    = get_field(av, 'classification',    'indeterminate');
    r.sizeCounts        = get_field(av, 'sizeCounts',        zeros(1, 0));
    r.lifetimeCounts    = get_field(av, 'lifetimeCounts',    zeros(1, 0));
    r.sethnaResidual    = get_field(av, 'sethnaResidual',    NaN);
    r.binMs             = get_field(av, 'binMs',             NaN);
    r.meanIeiMs         = get_field(av, 'meanIeiMs',         NaN);
    r.nAvalanches       = get_field(av, 'nAvalanches',       0);
    r.tailMassFrac      = get_field(av, 'tailMassFrac',      0);
end

function av = empty_av_struct()
    av = struct( ...
        'alphaLS',           NaN, ...
        'betaLS',            NaN, ...
        'alphaMLE',          NaN, ...
        'alphaMLE_ksp',      NaN, ...
        'lrtPreferredModel', 'power_law', ...
        'classification',    'indeterminate', ...
        'sizeCounts',        zeros(1, 0), ...
        'lifetimeCounts',    zeros(1, 0), ...
        'sethnaResidual',    NaN, ...
        'binMs',             NaN, ...
        'meanIeiMs',         NaN, ...
        'nAvalanches',       0, ...
        'tailMassFrac',      0);
end

function v = get_field(s, name, fallback)
    if isstruct(s) && isfield(s, name)
        v = s.(name);
    else
        v = fallback;
    end
end
