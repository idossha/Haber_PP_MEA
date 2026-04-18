function result = connectivity_xcorr(spikeTimesByChannel, varargin)
%CONNECTIVITY_XCORR Pairwise cross-correlation of binned spike trains.
%
%   result = CONNECTIVITY_XCORR(spikeTimesByChannel) bins each channel's
%   spike train at binMs (default 1 ms) and computes pairwise cross-
%   correlations within a [-maxLagMs, +maxLagMs] window. The peak
%   correlation across that window is taken as the connection weight; the
%   lag at which the peak occurs is also returned.
%
%   This is the cleaned-up successor to the ad-hoc cross-correlation logic
%   scattered across src/other_scripts/Ido_Ilhan__Correlation_automatic*.m
%   and OnlyCross_data_adaptive_zero_delay.m.
%
% INPUTS:
%   spikeTimesByChannel  -  Cell array (length nCh) of column vectors of
%                           spike timestamps in seconds. Empty cells become
%                           all-zero rows in the binned matrix.
%
% Name-value options:
%   'durationSec'    -  recording duration; default = max spike time across
%                       channels (rounded up to next bin).
%   'binMs'          -  bin width in ms (default 1).
%   'maxLagMs'       -  half-window of lags considered (default 100).
%   'normalization'  -  'pearson' (default) or 'coeff' (xcorr 'coeff') or
%                       'none'. 'pearson' rescales binned counts to z-scores
%                       so the peak is bounded in [-1, 1].
%
% OUTPUTS:
%   result  -  Struct with fields:
%       .adjacency      - nCh x nCh matrix of peak correlations
%       .peakLagMs      - nCh x nCh matrix of lags (ms) at the peak
%       .channels       - 1:nCh
%       .binMs          - bin width used
%       .maxLagMs       - lag half-window used
%       .normalization  - normalisation mode
%       .durationSec    - recording duration used
%
% Notes:
%   - The diagonal of adjacency is set to NaN (self-correlation excluded).
%   - The matrix is symmetric in magnitude (peakLag may differ in sign for
%     (i,j) vs (j,i)). Both halves are filled.

    p = inputParser;
    addRequired(p,  'spikeTimesByChannel', @iscell);
    addParameter(p, 'durationSec',    [], @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'binMs',          1,   @(x) isscalar(x) && x > 0);
    addParameter(p, 'maxLagMs',       100, @(x) isscalar(x) && x > 0);
    addParameter(p, 'normalization',  'pearson', ...
        @(s) any(strcmpi(s, {'pearson','coeff','none'})));
    parse(p, spikeTimesByChannel, varargin{:});
    opt = p.Results;

    nCh   = numel(spikeTimesByChannel);
    binSec = opt.binMs / 1000;

    % Recording duration: explicit value or inferred from spike data.
    if isempty(opt.durationSec)
        maxT = 0;
        for c = 1:nCh
            ts = spikeTimesByChannel{c};
            if ~isempty(ts)
                maxT = max(maxT, max(ts));
            end
        end
        durationSec = max(maxT, binSec);
    else
        durationSec = opt.durationSec;
    end
    nBins = max(1, floor(durationSec / binSec));

    % Build binned spike-count matrix (nCh x nBins).
    binned = zeros(nCh, nBins);
    edges  = (0:nBins) * binSec;
    for c = 1:nCh
        ts = spikeTimesByChannel{c};
        if isempty(ts); continue; end
        h = histcounts(ts, edges);
        binned(c, :) = h;
    end

    % Optional normalisation to make peaks comparable across channels.
    switch lower(opt.normalization)
        case 'pearson'
            mu = mean(binned, 2);
            sd = std(binned, 0, 2);
            sd(sd == 0) = 1;
            binnedZ = (binned - mu) ./ sd;
            xcorrMode = 'unbiased';
            scale = 1;  % xcorr 'unbiased' already normalizes by (N-|lag|)
        case 'coeff'
            binnedZ = binned;
            xcorrMode = 'coeff';
            scale = 1;
        case 'none'
            binnedZ = binned;
            xcorrMode = 'none';
            scale = 1;
    end

    maxLag = round(opt.maxLagMs / opt.binMs);
    adjacency = nan(nCh, nCh);
    peakLagMs = nan(nCh, nCh);

    for i = 1:nCh
        for j = i+1:nCh
            [c, lags] = xcorr(binnedZ(i, :), binnedZ(j, :), maxLag, xcorrMode);
            c = c * scale;
            [peakVal, idx] = max(c);
            adjacency(i, j) = peakVal;
            adjacency(j, i) = peakVal;
            peakLagMs(i, j) =  lags(idx) * opt.binMs;
            peakLagMs(j, i) = -lags(idx) * opt.binMs;
        end
    end

    result.adjacency     = adjacency;
    result.peakLagMs     = peakLagMs;
    result.channels      = 1:nCh;
    result.binMs         = opt.binMs;
    result.maxLagMs      = opt.maxLagMs;
    result.normalization = lower(opt.normalization);
    result.durationSec   = durationSec;
end
