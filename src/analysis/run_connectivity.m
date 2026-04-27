function results = run_connectivity(study, varargin)
%RUN_CONNECTIVITY Pairwise connectivity analysis across all pairs of a study.
%
%   results = RUN_CONNECTIVITY(study) loads spike times from cache for every
%   baseline/treatment pair in STUDY ('doi' or 'ket'), runs
%   connectivity_xcorr on each half, computes network_metrics on the
%   resulting adjacency matrices, and returns a struct array with one entry
%   per pair. No preprocessing is re-run; the cache must already exist
%   (build it with pipeline/preprocess_and_save.m).
%
%   RUN_CONNECTIVITY(study, 'save', true) also writes one analysis cache
%   file per pair into <cfg.paths.cache>/connectivity/ for fast re-loading.
%
% INPUTS:
%   study  -  'doi' or 'ket'.
%
% Name-value options:
%   'channels'       -  Channel list, default cfg.channels.default.
%   'binMs'          -  Bin width in ms, default 1.
%   'maxLagMs'       -  Lag half-window in ms, default 100.
%   'normalization'  -  'pearson' (default) | 'coeff' | 'none'.
%   'edgeThreshold'  -  Threshold for binarising the adjacency matrix in
%                       network_metrics, default 0.1.
%   'save'           -  Save per-pair connectivity caches, default false.
%
% OUTPUTS:
%   results  -  Struct array of length numel(pairs) with fields:
%       .pair             struct with .baseline, .treatment dataset names
%       .baseline         struct with .adjacency, .peakLagMs, .metrics
%       .treatment        struct with .adjacency, .peakLagMs, .metrics
%       .options          options used for connectivity_xcorr
%
% See also: CONNECTIVITY_XCORR, NETWORK_METRICS, LOAD_PAIR_SPIKES,
%           PREPROCESS_AND_SAVE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',      cfg.channels.default);
    addParameter(p, 'binMs',         cfg.connectivity.bin_ms,         @(x) isscalar(x) && x > 0);
    addParameter(p, 'maxLagMs',      cfg.connectivity.max_lag_ms,     @(x) isscalar(x) && x > 0);
    addParameter(p, 'normalization', cfg.connectivity.normalization,  ...
        @(s) any(strcmpi(s, {'pearson','coeff','none'})));
    % Default: target edge density via data-driven percentile thresholding.
    % An absolute threshold is rarely right because z-scored xcorr peak
    % magnitudes depend on bin size and recording duration.
    addParameter(p, 'edgeThreshold',  cfg.connectivity.edge_threshold, @(x) isempty(x) || isscalar(x));
    addParameter(p, 'edgeDensity',    cfg.connectivity.edge_density,   @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'save',           false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    pairs = get_pairs_and_labels(cfg, study);
    nPairs = numel(pairs);

    % Prepare save directory unconditionally so parfor can see it.
    cacheDir = fullfile(cfg.paths.cache, 'connectivity');
    if opt.save && ~exist(cacheDir, 'dir'); mkdir(cacheDir); end

    doSave   = opt.save;
    channels = opt.channels;
    binMs    = opt.binMs;
    maxLagMs = opt.maxLagMs;
    normMode = opt.normalization;
    edgeTh   = opt.edgeThreshold;
    edgeDen  = opt.edgeDensity;

    results = repmat(struct( ...
        'pair',      [], ...
        'baseline',  [], ...
        'treatment', [], ...
        'options',   []), 1, nPairs);

    % parfor runs in parallel when the Parallel Computing Toolbox is
    % available; otherwise it degrades gracefully to a serial for loop.
    parfor k = 1:nPairs
        [bSpikes, tSpikes, bMeta, tMeta] = ...
            load_pair_spikes(pairs(k), channels, cfg); %#ok<PFBNS>

        fprintf('run_connectivity(%s) [%d/%d] %s -> %s\n', ...
            study, k, nPairs, bMeta.datasetName, tMeta.datasetName);

        xcorrOpts = { ...
            'binMs',         binMs, ...
            'maxLagMs',      maxLagMs, ...
            'normalization', normMode};

        bRes = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, xcorrOpts{:});
        tRes = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, xcorrOpts{:});

        optK = struct('edgeThreshold', edgeTh, 'edgeDensity', edgeDen);
        bThresh = resolve_threshold(bRes.adjacency, optK);
        tThresh = resolve_threshold(tRes.adjacency, optK);

        bMetrics = network_metrics(bRes.adjacency, 'threshold', bThresh);
        tMetrics = network_metrics(tRes.adjacency, 'threshold', tThresh);

        r = struct( ...
            'pair',      pairs(k), ...
            'baseline',  struct( ...
                'datasetName', bMeta.datasetName, ...
                'adjacency',   bRes.adjacency, ...
                'peakLagMs',   bRes.peakLagMs, ...
                'metrics',     bMetrics), ...
            'treatment', struct( ...
                'datasetName', tMeta.datasetName, ...
                'adjacency',   tRes.adjacency, ...
                'peakLagMs',   tRes.peakLagMs, ...
                'metrics',     tMetrics), ...
            'options',   struct( ...
                'channels',           channels, ...
                'binMs',              binMs, ...
                'maxLagMs',           maxLagMs, ...
                'normalization',      normMode, ...
                'edgeThreshold',      edgeTh, ...
                'edgeDensity',        edgeDen, ...
                'thresholdBaseline',  bThresh, ...
                'thresholdTreatment', tThresh));
        results(k) = r;

        if doSave
            entry = r; %#ok<NASGU>
            outFile = fullfile(cacheDir, sprintf('connectivity_%s_%s.mat', ...
                study, sanitize(bMeta.datasetName)));
            save(outFile, '-struct', 'entry', '-v7.3');
            fprintf('           saved: %s\n', outFile);
        end
    end
end

% -------------------------------------------------------------------------
function s = sanitize(name)
    s = strrep(strrep(name, '-', '_'), '#', '_');
end

% -------------------------------------------------------------------------
function t = resolve_threshold(adjacency, opt)
% Explicit threshold wins. Otherwise pick the (1 - edgeDensity)
% percentile of the upper-triangular off-diagonal entries.
    if ~isempty(opt.edgeThreshold)
        t = opt.edgeThreshold;
        return;
    end
    n = size(adjacency, 1);
    upperMask = triu(true(n), 1);
    vals = adjacency(upperMask);
    vals = vals(~isnan(vals));
    if isempty(vals)
        t = 0;
        return;
    end
    pct = 100 * (1 - opt.edgeDensity);
    t = prctile(vals, pct);
end
