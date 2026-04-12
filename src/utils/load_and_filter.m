function [x, fs, durationSec, rawData] = load_and_filter(blockPath, channel, cfg)
%LOAD_AND_FILTER Read one TDT channel and apply the canonical filter stack.
%
%   [x, fs, durationSec, rawData] = LOAD_AND_FILTER(blockPath, channel, cfg)
%   loads channel CHANNEL from the TDT block at BLOCKPATH using TDTbin2mat,
%   then applies a 60 Hz notch stack (zero-phase, harmonics up to
%   cfg.filter.notch_max, Q = cfg.filter.notch_Q) followed by a zero-phase
%   Butterworth bandpass (cfg.filter.bp_low - cfg.filter.bp_high,
%   cfg.filter.bp_order). The TDT structure is also returned (with the
%   filtered data already substituted into the requested STORE), ready to
%   be passed straight to TDTthresh.
%
% INPUTS:
%   blockPath  -  Absolute path to a TDT block directory.
%   channel    -  Scalar channel index (1-based).
%   cfg        -  Project config struct from project_config().
%
% OUTPUTS:
%   x            -  Column vector of filtered samples (volts).
%   fs           -  Sampling rate (Hz) reported by TDT for this stream.
%   durationSec  -  Recording duration in seconds (numel(x) / fs).
%   rawData      -  Modified TDT struct: rawData.streams.(STORE).data has
%                   been replaced by the filtered signal as a row vector.
%                   This is the form expected by TDTthresh.
%
% Errors:
%   Re-throws TDTbin2mat errors so callers can decide whether to skip a
%   channel; the canonical preprocess pipeline catches and warns.

    store = cfg.tdt.store;

    rawData = TDTbin2mat(blockPath, 'STORE', store, 'CHANNEL', channel);
    fs = rawData.streams.(store).fs;
    x  = double(rawData.streams.(store).data(:));
    durationSec = numel(x) / fs;

    % Notch stack (zero-phase) at 60 Hz harmonics up to notch_max.
    notchStack = build_notch_filter(fs, ...
        cfg.filter.notch_f0, cfg.filter.notch_max, cfg.filter.notch_Q);
    for k = 1:numel(notchStack)
        x = filtfilt(notchStack(k).B, notchStack(k).A, x);
    end

    % Bandpass for spike band (zero-phase Butterworth).
    [bp_b, bp_a] = build_bandpass_filter(fs, ...
        cfg.filter.bp_low, cfg.filter.bp_high, cfg.filter.bp_order);
    x = filtfilt(bp_b, bp_a, x);

    % TDTthresh expects data as a row vector inside the TDT struct.
    rawData.streams.(store).data = x.';
end
