function verify_cache()
%VERIFY_CACHE Quick sanity summary of the v2.0 caches vs the literature.
%
%   Loads every cache file, prints per-dataset active-channel count
%   (baseline MFR >= cfg.silent.min_rate_spike), median spike rate on
%   active channels, and burst counts. Flags any dataset where spike
%   rates fall outside the Mossink 2021 / Brofiga 2023 reported ranges.
%
%   Run:
%     matlab -nodisplay -nosplash -batch \
%         "addpath(genpath('src'));addpath('scripts');verify_cache"

    cfg = project_config();

    allDs = [cfg.datasets.doi.baseline, cfg.datasets.doi.treatment, ...
             cfg.datasets.ket.baseline, cfg.datasets.ket.treatment];

    fprintf('\n==== verify_cache ====\n');
    fprintf('Literature-expected spike rate in rat cortical MEA at DIV ~7-21:\n');
    fprintf('   Mossink 2021  : median 0.1-5 spikes/s (~6-300 /min)\n');
    fprintf('   Chiappalone 2006: 0.5-10 /s (~30-600 /min)\n');
    fprintf('   Frega 2012    : 0.1-5 /s\n\n');
    fprintf('%-40s %6s %6s %8s %8s %6s\n', ...
        'dataset', 'fs', 'dur', 'actCh', 'medRate', 'nBrst');
    fprintf('%s\n', repmat('-', 1, 80));

    anomalies = {};
    for k = 1:numel(allDs)
        name = allDs{k};
        try
            c = load_cache(name, cfg);
        catch ME
            fprintf('%-40s  SKIP (%s)\n', name, ME.message);
            continue;
        end
        nActive = sum(c.spikeRates >= cfg.silent.min_rate_spike);
        medRate = median(c.spikeRates(c.spikeRates >= cfg.silent.min_rate_spike));
        nBurstTotal = sum(c.burstCounts);
        fprintf('%-40s %6.0f %6.0f %8d %8.2f %6d\n', ...
            name, c.fs, c.durationSec, nActive, medRate, nBurstTotal);

        if nActive < 5
            anomalies{end+1} = sprintf('%s: only %d active channels', name, nActive); %#ok<AGROW>
        end
        if medRate > 2000   % >33 Hz median is rare for primary cortical at DIV 7
            anomalies{end+1} = sprintf('%s: median rate %.0f /min is very high', name, medRate); %#ok<AGROW>
        end
    end

    fprintf('%s\n', repmat('-', 1, 80));
    if isempty(anomalies)
        fprintf('No anomalies flagged.\n');
    else
        fprintf('\n%d anomaly(ies):\n', numel(anomalies));
        for k = 1:numel(anomalies)
            fprintf('  - %s\n', anomalies{k});
        end
    end
    fprintf('\n');
end
