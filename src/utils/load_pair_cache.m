function [baselineCache, treatmentCache] = load_pair_cache(pair, cfg)
%LOAD_PAIR_CACHE Load both halves of a baseline/treatment pair.
%
%   [baselineCache, treatmentCache] = LOAD_PAIR_CACHE(pair, cfg) loads the
%   baseline and treatment cache structs for a single pair produced by
%   pair_datasets(). Both structs are full v2.0 caches (see LOAD_CACHE).
%
% INPUTS:
%   pair  -  Single element of the struct array from pair_datasets()
%            with fields .baseline and .treatment (dataset folder names).
%   cfg   -  Project config struct from project_config().
%
% OUTPUTS:
%   baselineCache   -  Cache struct for pair.baseline.
%   treatmentCache  -  Cache struct for pair.treatment.
%
% See also: LOAD_CACHE, LOAD_PAIR_METRIC, LOAD_PAIR_SPIKES.

    if nargin < 2
        error('load_pair_cache:Args', 'Usage: load_pair_cache(pair, cfg)');
    end
    if ~isfield(pair, 'baseline') || ~isfield(pair, 'treatment')
        error('load_pair_cache:Args', ...
            'pair must have .baseline and .treatment fields (see pair_datasets).');
    end

    baselineCache  = load_cache(pair.baseline,  cfg);
    treatmentCache = load_cache(pair.treatment, cfg);
end
