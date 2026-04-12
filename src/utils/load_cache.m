function cache = load_cache(datasetName, cfg)
%LOAD_CACHE Load a single preprocessed dataset cache with schema validation.
%
%   cache = LOAD_CACHE(datasetName, cfg) reads the cache file produced by
%   pipeline/preprocess_and_save.m for DATASETNAME. The returned struct
%   exposes the full v2.0 schema (see docs/CACHE_SCHEMA.md) with fields:
%
%     Metadata:
%       .version          char
%       .datasetName      char
%       .createdAt        char
%       .fs               double
%       .durationSec      double
%       .channelsUsed     (1 x nCh) int
%       .cfgUsed          struct
%
%     Spike-level:
%       .spikeTimes       (1 x nCh) cell of column vectors (seconds)
%       .spikeCounts      (nCh x 1) double
%       .spikeRates       (nCh x 1) double   (spikes/min)
%
%     Burst-level:
%       .burstStartIdx    (1 x nCh) cell, indices into spikeTimes{c}
%       .burstEndIdx      (1 x nCh) cell
%       .burstStartTimes  (1 x nCh) cell (seconds)
%       .burstEndTimes    (1 x nCh) cell (seconds)
%       .burstCounts      (nCh x 1) double
%       .burstRates       (nCh x 1) double   (bursts/min)
%
%   The loader performs:
%     1. File-exists check with an actionable error pointing at the
%        preprocessing script.
%     2. Schema-version sanity check: errors on unknown versions, warns
%        once on old versions.
%     3. Presence check on the canonical fields listed above.
%
% INPUTS:
%   datasetName  -  Dataset folder name (e.g. 'IdoDOI-230914-142502_#1').
%   cfg          -  Project config struct from project_config().
%
% OUTPUTS:
%   cache  -  Struct containing all v2.0 fields.
%
% Errors:
%   load_cache:MissingCache        - cache file not found on disk
%   load_cache:UnknownSchema       - cache version not recognised
%   load_cache:MissingField        - required v2.0 field missing from file
%
% See also: PREPROCESS_AND_SAVE, LOAD_PAIR_METRIC, LOAD_PAIR_SPIKES,
%           LOAD_PAIR_CACHE.

    if nargin < 2
        error('load_cache:Args', 'Usage: load_cache(datasetName, cfg)');
    end

    cachePath = fullfile(cfg.paths.cache, cache_filename(datasetName));
    if ~exist(cachePath, 'file')
        error('load_cache:MissingCache', ...
            ['No cache found for "%s".\n', ...
             'Expected: %s\n', ...
             'Run pipeline/preprocess_and_save.m first.'], ...
            datasetName, cachePath);
    end

    cache = load(cachePath);

    % Schema version dispatch.
    if ~isfield(cache, 'version')
        error('load_cache:UnknownSchema', ...
            ['Cache "%s" has no version field (likely legacy v1.x).\n', ...
             'Re-run pipeline/preprocess_and_save.m(''overwrite'', true) ', ...
             'to rebuild with the v2.0 schema.'], cachePath);
    end

    switch cache.version
        case '2.0'
            % OK - validate required fields below.
        otherwise
            error('load_cache:UnknownSchema', ...
                ['Cache "%s" has unknown schema version "%s".\n', ...
                 'This loader only understands v2.0. ', ...
                 'Re-run pipeline/preprocess_and_save.m(''overwrite'', true).'], ...
                cachePath, cache.version);
    end

    required = { ...
        'datasetName', 'createdAt', 'fs', 'durationSec', 'channelsUsed', ...
        'cfgUsed', ...
        'spikeTimes', 'spikeCounts', 'spikeRates', ...
        'burstStartIdx', 'burstEndIdx', 'burstStartTimes', 'burstEndTimes', ...
        'burstCounts', 'burstRates'};
    missing = setdiff(required, fieldnames(cache));
    if ~isempty(missing)
        error('load_cache:MissingField', ...
            ['Cache "%s" is missing required v2.0 field(s): %s\n', ...
             'Re-run pipeline/preprocess_and_save.m(''overwrite'', true).'], ...
            cachePath, strjoin(missing, ', '));
    end
end
