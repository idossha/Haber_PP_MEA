function [b, a] = build_bandpass_filter(fs, fLow, fHigh, order)
%BUILD_BANDPASS_FILTER Design a Butterworth bandpass for spike-band filtering.
%
%   [b, a] = BUILD_BANDPASS_FILTER(fs, fLow, fHigh, order) returns the
%   numerator and denominator coefficients of an order-N (default 4)
%   Butterworth bandpass with cutoffs [fLow, fHigh] in Hz, normalised to
%   the Nyquist frequency. Designed for use with FILTFILT (zero-phase),
%   matching the proof-of-concept pipeline (300-2500 Hz, order 4).
%
% INPUTS:
%   fs     -  Sampling rate in Hz.
%   fLow   -  Lower cutoff in Hz.
%   fHigh  -  Upper cutoff in Hz.
%   order  -  (optional) Filter order. Default = 4.
%
% OUTPUTS:
%   b, a   -  Butterworth bandpass coefficients.

    if nargin < 4 || isempty(order)
        order = 4;
    end
    if fs <= 0 || fLow <= 0 || fHigh <= 0 || fHigh <= fLow
        error('build_bandpass_filter:Args', ...
            'Require fs > 0 and 0 < fLow < fHigh.');
    end
    nyq = fs / 2;
    if fHigh >= nyq
        error('build_bandpass_filter:Args', ...
            'fHigh (%g Hz) must be < Nyquist (%g Hz).', fHigh, nyq);
    end
    [b, a] = butter(order, [fLow, fHigh] / nyq, 'bandpass');
end
