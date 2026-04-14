function result = connectivity_sttc(spikeTimesByChannel, varargin)
%CONNECTIVITY_STTC Pairwise spike-time tiling coefficient (STTC).
%
%   result = CONNECTIVITY_STTC(spikeTimesByChannel) computes the STTC
%   (Cutts & Eglen 2014) between every pair of channels. STTC is a
%   rate-independent correlation measure bounded in [-1, 1] that avoids the
%   systematic dependence on firing rate that affects binned Pearson
%   cross-correlation (Cutts & Eglen 2014, Table 3).
%
% INPUTS:
%   spikeTimesByChannel  -  Cell array (length nCh) of column vectors of
%                           spike timestamps in seconds. Empty cells are
%                           treated as silent channels.
%
% Name-value options:
%   'durationSec'  -  Total recording duration in seconds. Required for
%                     computing the fraction of time tiled by each train.
%                     Default = max spike time across channels.
%   'deltaMs'      -  Synchrony half-window in ms (default 50). Two spikes
%                     are considered synchronous if they fall within +/-
%                     deltaMs of each other. Cutts & Eglen (2014) recommend
%                     deltaMs matching the timescale of interest; 50 ms is
%                     a common choice for cortical MEA bursting dynamics.
%
% OUTPUTS:
%   result  -  Struct with fields:
%       .adjacency     - nCh x nCh symmetric matrix of STTC values
%       .channels      - 1:nCh
%       .deltaMs       - synchrony window used
%       .durationSec   - recording duration used
%
% REFERENCE:
%   Cutts CS, Eglen SJ (2014). Detecting pairwise correlations in spike
%   trains: an objective comparison of methods and application to the study
%   of retinal waves. J Neurosci 34:14288-14303.
%
% See also: CONNECTIVITY_XCORR, NETWORK_METRICS.

    p = inputParser;
    addRequired(p,  'spikeTimesByChannel', @iscell);
    addParameter(p, 'durationSec', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'deltaMs',    50,  @(x) isscalar(x) && x > 0);
    parse(p, spikeTimesByChannel, varargin{:});
    opt = p.Results;

    nCh     = numel(spikeTimesByChannel);
    deltaSec = opt.deltaMs / 1000;

    % Recording duration.
    if isempty(opt.durationSec)
        maxT = 0;
        for c = 1:nCh
            ts = spikeTimesByChannel{c};
            if ~isempty(ts)
                maxT = max(maxT, max(ts));
            end
        end
        durationSec = max(maxT, deltaSec);
    else
        durationSec = opt.durationSec;
    end

    % Pre-compute T_A for each channel: fraction of [0, durationSec] that
    % falls within +/- delta of any spike in the train.
    T = zeros(nCh, 1);
    for c = 1:nCh
        T(c) = tile_fraction(spikeTimesByChannel{c}, deltaSec, durationSec);
    end

    % Compute pairwise STTC.
    adjacency = nan(nCh, nCh);
    for i = 1:nCh
        spA = spikeTimesByChannel{i};
        for j = i+1:nCh
            spB = spikeTimesByChannel{j};
            sttc_val = compute_sttc(spA, spB, T(i), T(j), deltaSec);
            adjacency(i, j) = sttc_val;
            adjacency(j, i) = sttc_val;
        end
    end

    result.adjacency    = adjacency;
    result.channels     = 1:nCh;
    result.deltaMs      = opt.deltaMs;
    result.durationSec  = durationSec;
end

%% ========================================================================
function T = tile_fraction(spikeTimes, delta, duration)
%TILE_FRACTION Fraction of recording tiled by +/- delta windows around spikes.
%   Merges overlapping windows to avoid double-counting.
    if isempty(spikeTimes) || duration <= 0
        T = 0;
        return;
    end
    spikeTimes = sort(spikeTimes(:));
    nSp = numel(spikeTimes);

    % Build merged intervals [lo, hi] clipped to [0, duration].
    totalTiled = 0;
    lo = max(0, spikeTimes(1) - delta);
    hi = min(duration, spikeTimes(1) + delta);
    for k = 2:nSp
        newLo = max(0, spikeTimes(k) - delta);
        newHi = min(duration, spikeTimes(k) + delta);
        if newLo <= hi
            % Overlapping — extend current interval.
            hi = max(hi, newHi);
        else
            % Gap — close current interval and start new one.
            totalTiled = totalTiled + (hi - lo);
            lo = newLo;
            hi = newHi;
        end
    end
    totalTiled = totalTiled + (hi - lo);
    T = totalTiled / duration;
end

%% ========================================================================
function val = compute_sttc(spA, spB, TA, TB, delta)
%COMPUTE_STTC STTC between two spike trains (Cutts & Eglen 2014, eq. 1).
    nA = numel(spA);
    nB = numel(spB);

    % Handle degenerate cases.
    if nA == 0 || nB == 0
        val = 0;
        return;
    end

    % P_A: proportion of spikes in A within +/- delta of any spike in B.
    PA = count_near(spA, spB, delta) / nA;
    % P_B: proportion of spikes in B within +/- delta of any spike in A.
    PB = count_near(spB, spA, delta) / nB;

    % STTC = 0.5 * [ (PA - TB)/(1 - PA*TB) + (PB - TA)/(1 - PB*TA) ]
    % Guard against edge case where denominator = 0 (T = 1 and P = 1).
    denom1 = 1 - PA * TB;
    denom2 = 1 - PB * TA;
    if abs(denom1) < eps && abs(denom2) < eps
        val = 1;
        return;
    end
    term1 = 0;
    term2 = 0;
    if abs(denom1) >= eps
        term1 = (PA - TB) / denom1;
    end
    if abs(denom2) >= eps
        term2 = (PB - TA) / denom2;
    end
    val = 0.5 * (term1 + term2);
end

%% ========================================================================
function n = count_near(spA, spB, delta)
%COUNT_NEAR Count spikes in A that fall within +/- delta of any spike in B.
%   Uses a sorted-merge approach for O(nA + nB) efficiency.
    spA = sort(spA(:));
    spB = sort(spB(:));
    nA = numel(spA);
    nB = numel(spB);
    n  = 0;
    jLo = 1;
    for i = 1:nA
        tA = spA(i);
        % Advance jLo past spikes in B that are too early.
        while jLo <= nB && spB(jLo) < tA - delta
            jLo = jLo + 1;
        end
        % Check if any spike in B from jLo onward is within [tA-delta, tA+delta].
        if jLo <= nB && spB(jLo) <= tA + delta
            n = n + 1;
        end
    end
end
