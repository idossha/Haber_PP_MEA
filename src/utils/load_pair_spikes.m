function [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, channels, cfg)
%LOAD_PAIR_SPIKES Load spike-time cell arrays for a baseline/treatment pair.
%
%   [bSpikes, tSpikes, bMeta, tMeta] = LOAD_PAIR_SPIKES(pair, channels, cfg)
%   returns channel-aligned cell arrays of spike timestamps (seconds) for
%   both halves of a pair, restricted to CHANNELS. This is the canonical
%   entry point for connectivity, raster, and ISI analyses that operate on
%   raw spike times.
%
%   Output cell arrays are length numel(channels). For each requested
%   channel c, bSpikes{c} is a column vector of spike times from the
%   baseline recording, and tSpikes{c} is the corresponding vector from
%   the treatment recording. Missing channels (not present in
%   channelsUsed of either cache) receive empty vectors.
%
% INPUTS:
%   pair      -  Single element of the struct array from pair_datasets()
%                with .baseline and .treatment dataset folder names.
%   channels  -  Channel index list, e.g. 1:64. Output ordering follows
%                this list.
%   cfg       -  Project config struct from project_config().
%
% OUTPUTS:
%   bSpikes  -  (1 x numel(channels)) cell of column vectors (seconds).
%   tSpikes  -  (1 x numel(channels)) cell of column vectors (seconds).
%   bMeta    -  Struct with baseline metadata:
%                 .datasetName, .fs, .durationSec, .channels (the
%                 requested channels), .channelFound (logical same length
%                 as channels)
%   tMeta    -  Same struct for treatment.
%
% See also: LOAD_CACHE, LOAD_PAIR_CACHE, CONNECTIVITY_XCORR.

    if nargin < 3
        error('load_pair_spikes:Args', ...
            'Usage: load_pair_spikes(pair, channels, cfg)');
    end
    channels = channels(:)';
    nCh      = numel(channels);

    [loadedB, loadedT] = load_pair_cache(pair, cfg);

    bSpikes = repmat({zeros(0, 1)}, 1, nCh);
    tSpikes = repmat({zeros(0, 1)}, 1, nCh);

    bFound = false(1, nCh);
    tFound = false(1, nCh);

    chUsedB = loadedB.channelsUsed(:)';
    chUsedT = loadedT.channelsUsed(:)';
    [~, idxB] = ismember(channels, chUsedB);
    [~, idxT] = ismember(channels, chUsedT);

    for k = 1:nCh
        if idxB(k) > 0
            v = loadedB.spikeTimes{idxB(k)};
            bSpikes{k} = v(:);
            bFound(k)  = true;
        end
        if idxT(k) > 0
            v = loadedT.spikeTimes{idxT(k)};
            tSpikes{k} = v(:);
            tFound(k)  = true;
        end
    end

    bMeta = struct( ...
        'datasetName',  loadedB.datasetName, ...
        'fs',           loadedB.fs, ...
        'durationSec',  loadedB.durationSec, ...
        'channels',     channels, ...
        'channelFound', bFound);
    tMeta = struct( ...
        'datasetName',  loadedT.datasetName, ...
        'fs',           loadedT.fs, ...
        'durationSec',  loadedT.durationSec, ...
        'channels',     channels, ...
        'channelFound', tFound);
end
