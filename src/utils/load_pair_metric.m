function [bAll, tAll, chAll] = load_pair_metric(pairs, channels, metric, cfg)
%LOAD_PAIR_METRIC Pool a per-channel summary metric across all pairs.
%
%   [bAll, tAll, chAll] = LOAD_PAIR_METRIC(pairs, channels, metric, cfg)
%   loads cached summary-metric vectors from each pair, restricts to the
%   requested CHANNELS, drops NaNs, and concatenates the channel
%   observations across pairs. The returned vectors are aligned: bAll(k),
%   tAll(k), and chAll(k) refer to the same channel observation.
%
% INPUTS:
%   pairs     -  Struct array from pair_datasets().
%   channels  -  Channel index list (e.g. 1:64).
%   metric    -  Name of the per-channel summary field to pool. One of:
%                   'spikeRates'   (spikes/min)
%                   'spikeCounts'  (total spike count)
%                   'burstRates'   (bursts/min)
%                   'burstCounts'  (total burst count)
%   cfg       -  Project config struct from project_config().
%
% OUTPUTS:
%   bAll   -  Pooled baseline values (column vector).
%   tAll   -  Pooled treatment values (column vector).
%   chAll  -  Channel index for each observation (row vector).
%
% Errors:
%   load_pair_metric:UnknownMetric - metric not in the allowed set
%   load_cache:MissingCache         - cache file missing on disk
%
% See also: LOAD_CACHE, LOAD_PAIR_CACHE, LOAD_PAIR_SPIKES.

    if nargin < 4
        error('load_pair_metric:Args', ...
            'Usage: load_pair_metric(pairs, channels, metric, cfg)');
    end
    allowed = {'spikeRates', 'spikeCounts', 'burstRates', 'burstCounts'};
    if ~ismember(metric, allowed)
        error('load_pair_metric:UnknownMetric', ...
            'metric must be one of {%s}, got "%s".', ...
            strjoin(allowed, ', '), metric);
    end

    bAll  = [];
    tAll  = [];
    chAll = [];

    for p = 1:numel(pairs)
        [loadedB, loadedT] = load_pair_cache(pairs(p), cfg);

        valuesB = loadedB.(metric)(:);
        valuesT = loadedT.(metric)(:);
        chUsedB = loadedB.channelsUsed(:)';
        chUsedT = loadedT.channelsUsed(:)';

        [~, idxB] = ismember(channels, chUsedB);
        [~, idxT] = ismember(channels, chUsedT);
        validIdx = (idxB > 0) & (idxT > 0);
        if ~any(validIdx)
            continue;
        end

        bPair  = valuesB(idxB(validIdx));
        tPair  = valuesT(idxT(validIdx));
        chPair = channels(validIdx);

        keep   = ~isnan(bPair) & ~isnan(tPair);
        bPair  = bPair(keep);
        tPair  = tPair(keep);
        chPair = chPair(keep);

        bAll  = [bAll;  bPair];   %#ok<AGROW>
        tAll  = [tAll;  tPair];   %#ok<AGROW>
        chAll = [chAll, chPair];  %#ok<AGROW>
    end
end
