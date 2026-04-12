function notchStack = build_notch_filter(fs, f0, fmax, Q)
%BUILD_NOTCH_FILTER Pre-compute IIR notch coefficients at f0 and harmonics.
%
%   notchStack = BUILD_NOTCH_FILTER(fs, f0, fmax, Q) returns a struct array
%   with one entry per notch (fundamental + integer harmonics) suitable for
%   sequential zero-phase filtering with FILTFILT. Each entry has fields:
%       .f      - centre frequency (Hz)
%       .B, .A  - IIR notch numerator/denominator coefficients (length 3)
%
%   The harmonics included are f0, 2*f0, ..., up to but not exceeding both
%   fmax and the Nyquist frequency (fs/2). This matches the canonical 60 Hz
%   stack used in the proof-of-concept figure scripts (f0 = 60, fmax = 780,
%   Q = 35).
%
% INPUTS:
%   fs    -  Sampling rate in Hz.
%   f0    -  Fundamental notch frequency in Hz (e.g. 60).
%   fmax  -  Highest harmonic to include in Hz (e.g. 780).
%   Q     -  Quality factor of each notch (higher = narrower).
%
% OUTPUTS:
%   notchStack  -  Struct array, length = number of notches actually built.
%
% Notes:
%   - This function intentionally has **no toolbox dependency**. The
%     coefficients come from the standard second-order IIR notch (the
%     Robert Bristow-Johnson "Audio EQ Cookbook" biquad notch, which is
%     the same recipe used internally by DSP System Toolbox's iirnotch
%     / designNotchPeakIIR). We implement it directly so the pipeline
%     runs with only Signal Processing Toolbox (needed elsewhere for
%     filtfilt, butter and xcorr).
%   - Applied zero-phase with filtfilt so the effective magnitude
%     response is squared; this is by design and matches the POC.

    if nargin < 4
        error('build_notch_filter:Args', ...
            'Usage: build_notch_filter(fs, f0, fmax, Q)');
    end
    if fs <= 0 || f0 <= 0 || fmax <= 0 || Q <= 0
        error('build_notch_filter:Args', ...
            'fs, f0, fmax, and Q must all be strictly positive.');
    end

    nyq   = fs / 2;
    freqs = f0:f0:fmax;
    freqs = freqs(freqs < nyq);

    notchStack = repmat(struct('f', [], 'B', [], 'A', []), 1, numel(freqs));
    for k = 1:numel(freqs)
        fk = freqs(k);
        [B, A] = biquad_notch(fk, fs, Q);
        notchStack(k).f = fk;
        notchStack(k).B = B;
        notchStack(k).A = A;
    end
end

% =========================================================================
function [B, A] = biquad_notch(f0, fs, Q)
% Standard second-order IIR notch (RBJ biquad cookbook):
%
%     w0 = 2*pi*f0/fs
%     alpha = sin(w0) / (2*Q)
%
%     b0 =  1         a0 =  1 + alpha
%     b1 = -2*cos(w0) a1 = -2*cos(w0)
%     b2 =  1         a2 =  1 - alpha
%
%   Normalised so a0 = 1.
%
% This is identical to the coefficients returned by iirnotch(w0/pi, bw/pi)
% when bw = w0/Q, and matches designNotchPeakIIR(..., QualityFactor=Q).

    w0    = 2 * pi * f0 / fs;
    alpha = sin(w0) / (2 * Q);
    cosw0 = cos(w0);

    b0 = 1;
    b1 = -2 * cosw0;
    b2 = 1;
    a0 = 1 + alpha;
    a1 = -2 * cosw0;
    a2 = 1 - alpha;

    B = [b0, b1, b2] / a0;
    A = [1,  a1/a0, a2/a0];
end
