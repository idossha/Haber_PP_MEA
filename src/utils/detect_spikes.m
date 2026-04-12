function spikeTimes = detect_spikes(filteredTdtStruct, cfg)
%DETECT_SPIKES Wrap TDTthresh with the canonical project parameters.
%
%   spikeTimes = DETECT_SPIKES(filteredTdtStruct, cfg) runs TDTthresh on a
%   pre-filtered TDT data structure (as returned by load_and_filter) and
%   returns the sorted spike timestamps in seconds. If TDTthresh returns no
%   detections (or any field is missing), an empty column vector is
%   returned.
%
% INPUTS:
%   filteredTdtStruct  -  TDT data struct with filtered samples already
%                         substituted into .streams.(cfg.tdt.store).data.
%   cfg                -  Project config struct from project_config().
%
% OUTPUTS:
%   spikeTimes  -  Sorted column vector of spike timestamps in seconds.
%
% Detection parameters (from cfg.spike): MODE, POLARITY, STD, TAU.

    store = cfg.tdt.store;

    % cfg.spike.tau is the TDTthresh noise-window length in SECONDS
    % (not a refractory period; see TDTthresh.m line 35-36). Passing it
    % verbatim restores the POC / TDT-default behaviour.
    spikes = TDTthresh(filteredTdtStruct, store, ...
        'MODE',     cfg.spike.mode, ...
        'POLARITY', cfg.spike.polarity, ...
        'STD',      cfg.spike.std_factor, ...
        'TAU',      cfg.spike.tau);

    if isempty(spikes) ...
            || ~isfield(spikes, 'snips') ...
            || ~isfield(spikes.snips, 'Snip') ...
            || ~isfield(spikes.snips.Snip, 'ts')
        spikeTimes = zeros(0, 1);
        return;
    end

    spikeTimes = sort(spikes.snips.Snip.ts(:));

    % Post-detection refractory enforcement. TDTthresh's TAU is a noise
    % window, not a refractory period; several literature references
    % (Brofiga 2023 p.4, Garofalo 2009 p.3) enforce a 1-2 ms refractory
    % window after detection. We apply it here by dropping any spike
    % within refractory_ms of the previous retained spike.
    refMs = 0;
    if isfield(cfg.spike, 'refractory_ms') && ~isempty(cfg.spike.refractory_ms)
        refMs = cfg.spike.refractory_ms;
    end
    if refMs > 0 && numel(spikeTimes) >= 2
        refSec = refMs * 1e-3;
        keep = true(size(spikeTimes));
        lastKept = spikeTimes(1);
        for k = 2:numel(spikeTimes)
            if spikeTimes(k) - lastKept < refSec
                keep(k) = false;
            else
                lastKept = spikeTimes(k);
            end
        end
        spikeTimes = spikeTimes(keep);
    end
end
