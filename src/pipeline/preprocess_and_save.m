function preprocess_and_save(varargin)
%PREPROCESS_AND_SAVE Build per-dataset spike/burst caches for all studies.
%
%   PREPROCESS_AND_SAVE() runs the canonical preprocessing pipeline for
%   every dataset listed in project_config(): load TDT block -> notch +
%   bandpass filter -> TDTthresh spike detection -> ISI burst detection ->
%   per-channel spike times, burst indices, rates, and counts -> save one
%   .mat file per dataset under cfg.paths.cache. This is the only place
%   that touches raw TDT data; every figure and analysis script downstream
%   loads from cache via utils/load_cache.m.
%
%   PREPROCESS_AND_SAVE('overwrite', true) re-runs even when a cache file
%   already exists. By default existing caches are skipped.
%
%   PREPROCESS_AND_SAVE('study', 'doi')   processes only the DOI study.
%   PREPROCESS_AND_SAVE('study', 'ket')   processes only the ketanserin study.
%   PREPROCESS_AND_SAVE('study', 'all')   (default) processes both.
%
% INPUTS:
%   Name-value options:
%     'overwrite'  -  logical, default false
%     'study'      -  one of {'all','doi','ket'}, default 'all'
%
% OUTPUTS:
%   One cache file per dataset at:
%       <cfg.paths.cache>/cache_<safeDatasetName>.mat
%
%   See docs/CACHE_SCHEMA.md for the full on-disk schema (v2.0). Summary:
%
%     Metadata:
%       version          char  ('2.0')
%       datasetName      char
%       createdAt        char  (ISO 8601)
%       fs               double
%       durationSec      double
%       channelsUsed     (1 x nCh) int
%       cfgUsed          struct (snapshot of project_config at build time)
%
%     Spike-level (new in v2.0):
%       spikeTimes       (1 x nCh) cell, each cell a column vector of
%                        timestamps in seconds (sorted ascending)
%
%     Spike summary:
%       spikeCounts      (nCh x 1) double   - total spikes per channel
%       spikeRates       (nCh x 1) double   - spikes/min per channel
%
%     Burst-level (new in v2.0):
%       burstStartIdx    (1 x nCh) cell, each cell a row vector of start
%                        indices into spikeTimes{c} (empty if no bursts)
%       burstEndIdx      (1 x nCh) cell, same shape, end indices
%       burstStartTimes  (1 x nCh) cell, row vector of start times (s)
%       burstEndTimes    (1 x nCh) cell, row vector of end times (s)
%
%     Burst summary:
%       burstCounts      (nCh x 1) double   - total bursts per channel
%       burstRates       (nCh x 1) double   - bursts/min per channel
%
% Notes:
%   - Filtering and spike detection parameters come from cfg.filter and
%     cfg.spike. Edit src/config/project_config.m to change them.
%   - Burst detection uses cfg.burst (ISI-threshold; matches POC
%     burst_stats.m verbatim).
%   - Rates are NaN if no channel produced a valid recording duration.
%   - On a per-channel load/detect failure, the channel's slots are left
%     NaN / empty and a warning is emitted; preprocessing continues.

    p = inputParser;
    addParameter(p, 'overwrite', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'study',     'all', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});
    overwrite = p.Results.overwrite;
    study     = lower(char(p.Results.study));

    cfg = project_config();
    addpath(genpath(cfg.paths.tdt_sdk));

    if ~exist(cfg.paths.cache, 'dir')
        mkdir(cfg.paths.cache);
    end

    % Build the flat list of datasets to process for the requested study.
    datasets = {};
    if any(strcmp(study, {'all','doi'}))
        datasets = [datasets, cfg.datasets.doi.baseline, cfg.datasets.doi.treatment];
    end
    if any(strcmp(study, {'all','ket'}))
        datasets = [datasets, cfg.datasets.ket.baseline, cfg.datasets.ket.treatment];
    end
    if isempty(datasets)
        error('preprocess_and_save:NoDatasets', ...
            'No datasets selected for study "%s".', study);
    end

    nTotal = numel(datasets);
    fprintf('Preprocessing %d datasets (study=%s, overwrite=%d)\n', ...
        nTotal, study, overwrite);

    for d = 1:nTotal
        datasetName = datasets{d};
        cachePath   = fullfile(cfg.paths.cache, cache_filename(datasetName));
        if exist(cachePath, 'file') && ~overwrite
            fprintf('  [%d/%d] cache exists, skipping: %s\n', d, nTotal, datasetName);
            continue;
        end

        try
            blockPath = dataset_path(datasetName, cfg);
        catch ME
            warning('preprocess_and_save:DatasetMissing', ...
                'Skipping "%s": %s', datasetName, ME.message);
            continue;
        end

        fprintf('  [%d/%d] processing: %s\n', d, nTotal, datasetName);
        cacheStruct = process_one_dataset(blockPath, cfg);

        cacheStruct.version     = '2.0';
        cacheStruct.datasetName = datasetName;
        cacheStruct.createdAt   = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
        cacheStruct.cfgUsed     = cfg;

        save(cachePath, '-struct', 'cacheStruct', '-v7.3');
        fprintf('         saved: %s\n', cachePath);
    end

    fprintf('Preprocessing complete. Cache directory: %s\n', cfg.paths.cache);
end

% =========================================================================
function result = process_one_dataset(blockPath, cfg)
% Process a single TDT block. For every channel listed in
% cfg.channels.default: load + filter -> detect spikes -> detect bursts ->
% compute rates. Returns the cache struct (without metadata fields).

    channels = cfg.channels.default;
    nCh      = numel(channels);

    % Per-channel storage.
    spikeTimes      = cell(1, nCh);
    spikeCounts     = NaN(nCh, 1);
    spikeRates      = NaN(nCh, 1);
    burstStartIdx   = cell(1, nCh);
    burstEndIdx     = cell(1, nCh);
    burstStartTimes = cell(1, nCh);
    burstEndTimes   = cell(1, nCh);
    burstCounts     = NaN(nCh, 1);
    burstRates      = NaN(nCh, 1);

    fs          = NaN;
    durationSec = NaN;

    wb = waitbar(0, 'Processing channels...');
    cleanupWb = onCleanup(@() safeClose(wb));

    for i = 1:nCh
        ch = channels(i);
        try
            [~, fsCh, durCh, filteredStruct] = load_and_filter(blockPath, ch, cfg);
        catch ME
            warning('preprocess_and_save:ChannelLoad', ...
                'Channel %d failed to load: %s', ch, ME.message);
            spikeTimes{i}      = zeros(0, 1);
            burstStartIdx{i}   = zeros(1, 0);
            burstEndIdx{i}     = zeros(1, 0);
            burstStartTimes{i} = zeros(1, 0);
            burstEndTimes{i}   = zeros(1, 0);
            updateWb(wb, i, nCh, ch);
            continue;
        end

        if isnan(fs)
            fs          = fsCh;
            durationSec = durCh;
        end

        % ---- Spike detection ----
        ts = detect_spikes(filteredStruct, cfg);
        spikeTimes{i}  = ts;
        spikeCounts(i) = numel(ts);
        if durationSec > 0
            spikeRates(i) = spikeCounts(i) * (60 / durationSec);
        end

        % ---- Burst detection ----
        if numel(ts) >= cfg.burst.min_spikes
            [bStartIdx, bEndIdx, bCount] = detect_bursts(ts, cfg);
        else
            bStartIdx = zeros(1, 0);
            bEndIdx   = zeros(1, 0);
            bCount    = 0;
        end
        burstStartIdx{i} = bStartIdx;
        burstEndIdx{i}   = bEndIdx;

        if bCount > 0
            burstStartTimes{i} = ts(bStartIdx).';
            burstEndTimes{i}   = ts(bEndIdx).';
        else
            burstStartTimes{i} = zeros(1, 0);
            burstEndTimes{i}   = zeros(1, 0);
        end
        burstCounts(i) = bCount;
        if durationSec > 0
            burstRates(i) = (bCount / durationSec) * 60;
        end

        updateWb(wb, i, nCh, ch);
    end

    if isnan(durationSec) || durationSec <= 0
        % If no channel produced a valid duration, all per-minute rates
        % are undefined.
        spikeRates(:) = NaN;
        burstRates(:) = NaN;
    end

    result.fs              = fs;
    result.durationSec     = durationSec;
    result.channelsUsed    = channels;

    result.spikeTimes      = spikeTimes;
    result.spikeCounts     = spikeCounts;
    result.spikeRates      = spikeRates;

    result.burstStartIdx   = burstStartIdx;
    result.burstEndIdx     = burstEndIdx;
    result.burstStartTimes = burstStartTimes;
    result.burstEndTimes   = burstEndTimes;
    result.burstCounts     = burstCounts;
    result.burstRates      = burstRates;
end

% -------------------------------------------------------------------------
function updateWb(wb, i, n, ch)
    if isvalid(wb)
        waitbar(i / n, wb, sprintf('Processing channel %d (%d / %d)', ch, i, n));
    end
end

% -------------------------------------------------------------------------
function safeClose(wb)
    if isvalid(wb)
        close(wb);
    end
end
