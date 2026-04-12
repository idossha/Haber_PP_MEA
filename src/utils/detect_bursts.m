function [burstStarts, burstEnds, burstCount] = detect_bursts(spikeTimes, cfg)
%DETECT_BURSTS ISI-threshold burst detector (matches POC burst_stats.m).
%
%   [burstStarts, burstEnds, burstCount] = DETECT_BURSTS(spikeTimes, cfg)
%   labels bursts in a sorted spike-time vector using the same ISI state
%   machine as src/oldpipeline_scripts/burst_stats.m:
%
%       - A burst starts when ISI <= isi_max_ms.
%       - A burst ends when ISI >  isi_end_ms.
%       - A burst is kept only if it contains at least min_spikes spikes.
%       - Trailing bursts that run to the end of the recording are kept if
%         they meet the min_spikes rule.
%
%   Thresholds come from cfg.burst (isi_max_ms, isi_end_ms, min_spikes).
%
% INPUTS:
%   spikeTimes  -  Sorted column vector of spike timestamps (seconds).
%   cfg         -  Project config struct from project_config().
%
% OUTPUTS:
%   burstStarts  -  Row vector of start indices into spikeTimes for each
%                   detected burst (one entry per burst).
%   burstEnds    -  Row vector of end indices (same length as burstStarts).
%   burstCount   -  numel(burstStarts).

    isiMax    = cfg.burst.isi_max_ms / 1000;     % seconds
    isiEnd    = cfg.burst.isi_end_ms / 1000;     % seconds
    minSpikes = cfg.burst.min_spikes;

    burstStarts = zeros(1, 0);
    burstEnds   = zeros(1, 0);

    if numel(spikeTimes) < minSpikes
        burstCount = 0;
        return;
    end

    ISI     = diff(spikeTimes);
    inBurst = false;
    startIdx = NaN;

    for ii = 1:numel(ISI)
        if ~inBurst && ISI(ii) <= isiMax
            inBurst  = true;
            startIdx = ii;
        elseif inBurst && ISI(ii) > isiEnd
            endIdx = ii;
            if endIdx - startIdx + 1 >= minSpikes
                burstStarts(end+1) = startIdx; %#ok<AGROW>
                burstEnds(end+1)   = endIdx;   %#ok<AGROW>
            end
            inBurst = false;
        end
    end

    % Tail burst: if we are still inside a burst at the end of the trace,
    % close it at the last spike.
    if inBurst
        endIdx = numel(spikeTimes);
        if endIdx - startIdx + 1 >= minSpikes
            burstStarts(end+1) = startIdx;
            burstEnds(end+1)   = endIdx;
        end
    end

    burstCount = numel(burstStarts);
end
