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

    % --- Vectorised FFT-based cross-correlation ---
    % One FFT per channel (60 total) instead of two per pair (3540 total).
    % Cross-spectra are computed via element-wise multiply in frequency
    % domain, then batch-IFFT extracts only the ±maxLag window.
    % We only need circular xcorr to match linear xcorr at |lag| <= maxLag,
    % so nfft >= nBins + maxLag suffices (much smaller than 2*nBins - 1).
    nfft = 2^nextpow2(nBins + maxLag);
    F    = fft(binnedZ, nfft, 2);        % nCh x nfft  (all channels at once)

    lagRange = -maxLag : maxLag;         % e.g. -100:100
    nLags    = numel(lagRange);

    % IFFT lag-index mapping: lag 0 -> col 1, lag +k -> col k+1,
    % lag -k -> col nfft-k+1.
    posIdx = 1 : maxLag + 1;                      % lags 0..+maxLag
    negIdx = nfft - maxLag + 1 : nfft;            % lags -maxLag..-1
    extractIdx = [negIdx, posIdx];                 % ordered -maxLag..+maxLag

    % Normalisation denominator for each lag.
    switch lower(opt.normalization)
        case 'pearson'
            denom = nBins - abs(lagRange);         % 'unbiased': N-|k|
        case 'coeff'
            % 'coeff' divides by sqrt(Rxx(0)*Ryy(0)); handled per-pair below.
            denom = ones(1, nLags);
        otherwise
            denom = ones(1, nLags);
    end

    for i = 1:nCh - 1
        jIdx = (i + 1) : nCh;
        nJ   = numel(jIdx);

        % Cross-spectral density: conj(F_i) .* F_j  for all j > i.
        S = bsxfun(@times, conj(F(i, :)), F(jIdx, :));   % nJ x nfft

        % Batch IFFT and extract the ±maxLag window.
        cFull    = real(ifft(S, [], 2));                   % nJ x nfft
        cTrimmed = cFull(:, extractIdx);                   % nJ x nLags

        % Apply normalisation.
        if strcmpi(opt.normalization, 'coeff')
            % 'coeff': divide by sqrt(Rxx(0)*Ryy(0)) = sqrt(energy_i * energy_j)
            Ei  = sum(binnedZ(i, :).^2);
            Ej  = sum(binnedZ(jIdx, :).^2, 2);
            cTrimmed = cTrimmed ./ sqrt(Ei .* Ej);
        else
            cTrimmed = cTrimmed ./ denom;
        end

        % Peak value and lag for each pair.
        [peakVals, peakPos] = max(cTrimmed, [], 2);

        adjacency(i, jIdx) = peakVals;
        adjacency(jIdx, i) = peakVals;
        peakLagMs(i, jIdx) =  lagRange(peakPos) * opt.binMs;
        peakLagMs(jIdx, i) = -lagRange(peakPos) * opt.binMs;
    end

    result.adjacency     = adjacency;
    result.peakLagMs     = peakLagMs;
    result.channels      = 1:nCh;
    result.binMs         = opt.binMs;
    result.maxLagMs      = opt.maxLagMs;
    result.normalization = lower(opt.normalization);
    result.durationSec   = durationSec;
end
