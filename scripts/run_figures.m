function run_figures(varargin)
%RUN_FIGURES Single entry point for all figure generation.
%
%   run_figures()                                  All figures, both studies, all pairs.
%   run_figures('study', 'doi')                    DOI study only.
%   run_figures('study', 'ket')                    Ketanserin study only.
%   run_figures('category', 'rates')               Rate panels only.
%   run_figures('category', 'connectivity')        Connectivity panels only.
%   run_figures('category', 'si')                  Supplementary + shared panels.
%   run_figures('study','doi','category','rates')  DOI rates only.
%
%   Categories:
%     'rates'        — traces, rasters, scatter, violin, bar, bootstrap
%     'connectivity' — exemplar + summary connectivity panels
%     'si'           — supplementary panels, electrode map, correlograms
%
%   By default, per-pair panels (traces, rasters, connectivity) are
%   generated for EVERY baseline/treatment pair.  Pass 'allPairs', false
%   to generate only the default representative pair (faster, for testing).
%
%   Preconditions: pipeline/preprocess_and_save.m has been run so the
%   v2.0 caches exist under <repo>/cache/.
%
%   Output: 600 DPI PDF + PNG files under <repo>/output/.
%
% See also: FIG_CONNECTIVITY_EXEMPLAR, FIG_CONNECTIVITY_SUMMARY,
%           PROJECT_CONFIG.

    p = inputParser;
    addParameter(p, 'study',    'all', ...
        @(s) any(strcmpi(s, {'all','doi','ket'})));
    addParameter(p, 'category', 'all', ...
        @(s) any(strcmpi(s, {'all','rates','connectivity','si'})));
    addParameter(p, 'allPairs', true, ...
        @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    switch lower(opt.study)
        case 'all'; studies = {'doi','ket'};
        case 'doi'; studies = {'doi'};
        case 'ket'; studies = {'ket'};
    end

    cat            = lower(opt.category);
    doRates        = any(strcmp(cat, {'all','rates'}));
    doConnectivity = any(strcmp(cat, {'all','connectivity'}));
    doSI           = any(strcmp(cat, {'all','si'}));

    t0 = tic;
    nTotal = 0; nOk = 0; nErr = 0;
    errors = {};

    fprintf('\n==== run_figures ====\n');
    fprintf('  studies  : %s\n', strjoin(studies, ', '));
    fprintf('  category : %s\n', cat);
    fprintf('  allPairs : %d\n\n', opt.allPairs);

    cfg = project_config();

    % --- Per-study figures ------------------------------------------------
    for s = 1:numel(studies)
        study = studies{s};
        fprintf('--- %s ---\n', upper(study));

        % — Rate panels ----------------------------------------------------
        if doRates
            % Aggregate panels (pool all pairs internally)
            calls = {
                'fig_rate_scatter',               {'spike', study}
                'fig_rate_scatter',               {'burst', study}
                'fig_spike_rate',                 {study}
                'fig_burst_rate',                 {study}
                'fig_spike_pct_change_violin',    {study}
                'fig_burst_pct_change_violin',    {study}
                'fig_rate_change_bar',            {'spike', study}
                'fig_rate_change_bar',            {'burst', study}
                'fig_stats_bootstrap',            {study}
            };
            [nTotal, nOk, nErr, errors] = run_calls( ...
                calls, nTotal, nOk, nErr, errors);

            % Per-pair panels (traces, rasters)
            [pairs, ~] = get_pairs_and_labels(cfg, study);
            tracesDir = fullfile(output_path(cfg, study, 'rates', ''), 'traces');
            rasterDir = fullfile(output_path(cfg, study, 'rates', ''), 'raster');
            if opt.allPairs
                pairList = 1:numel(pairs);
            else
                pairList = max(1, round(numel(pairs) / 2));
            end
            for k = pairList
                tag = sprintf('%s_pair%d', study, k);
                calls = {
                    'fig_rate_panel_traces', ...
                        {study, 'pairIndex', k, ...
                         'outName', [tag '_traces'], 'outDir', tracesDir}
                    'fig_rate_panel_raster', ...
                        {study, 'pairIndex', k, ...
                         'outName', [tag '_raster'], 'outDir', rasterDir}
                };
                [nTotal, nOk, nErr, errors] = run_calls( ...
                    calls, nTotal, nOk, nErr, errors);
            end
        end

        % — Connectivity panels --------------------------------------------
        if doConnectivity
            [pairs, ~] = get_pairs_and_labels(cfg, study);
            if opt.allPairs
                pairList = 1:numel(pairs);
            else
                pairList = max(1, round(numel(pairs) / 2));
            end

            % Exemplar panels for each pair
            for k = pairList
                prefix = sprintf('%s_pair%d', study, k);
                calls = {
                    'fig_connectivity_exemplar', ...
                        {study, 'pairIndex', k, 'prefix', prefix}
                };
                [nTotal, nOk, nErr, errors] = run_calls( ...
                    calls, nTotal, nOk, nErr, errors);
            end

            % Pre-compute connectivity once for CDF + summary panels
            nTotal = nTotal + 1;
            connLabel = sprintf('run_connectivity(%s)', study);
            results_conn = [];
            try
                fprintf('  [%2d] %s\n', nTotal, connLabel);
                results_conn = run_connectivity(study);
                nOk = nOk + 1;
            catch ME
                nErr = nErr + 1;
                fprintf('       ERROR: %s\n', ME.message);
                errors{end+1} = sprintf('%s: %s', connLabel, ME.message); %#ok<AGROW>
                if ~isempty(ME.stack)
                    fprintf('       at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
                end
            end

            if ~isempty(results_conn)
                calls = {
                    'fig_connectivity_cdf',     {study, 'results', results_conn}
                    'fig_connectivity_summary', {study, 'results', results_conn}
                };
                [nTotal, nOk, nErr, errors] = run_calls( ...
                    calls, nTotal, nOk, nErr, errors);
            end

            % Python heatmaps from exported CSVs
            nTotal = nTotal + 1;
            label  = sprintf('plot_heatmaps.py --study %s', study);
            try
                fprintf('  [%2d] %s\n', nTotal, label);
                pyScript = fullfile(cfg.paths.root, 'scripts', 'plot_heatmaps.py');
                [status, result] = system(sprintf('python3 "%s" --study %s', pyScript, study));
                if status ~= 0
                    error('run_figures:PythonFailed', ...
                        'plot_heatmaps.py failed:\n%s', result);
                end
                fprintf('%s', result);
                nOk = nOk + 1;
            catch ME
                nErr = nErr + 1;
                fprintf('       ERROR: %s\n', ME.message);
                errors{end+1} = sprintf('%s: %s', label, ME.message); %#ok<AGROW>
            end
        end

        % — Supplementary panels -------------------------------------------
        if doSI
            calls = {
                'fig_si_zoomed_raster',          {study}
                'fig_si_burst_scatter',          {study}
                'fig_si_burst_overlay_raster',   {study}
            };
            [nTotal, nOk, nErr, errors] = run_calls( ...
                calls, nTotal, nOk, nErr, errors);

            % Python correlograms (CSV exported by fig_connectivity_exemplar)
            nTotal = nTotal + 1;
            label  = sprintf('plot_correlograms.py --study %s', study);
            try
                fprintf('  [%2d] %s\n', nTotal, label);
                pyScript = fullfile(cfg.paths.root, 'scripts', 'plot_correlograms.py');
                [status, result] = system(sprintf('python3 "%s" --study %s', pyScript, study));
                if status ~= 0
                    error('run_figures:PythonFailed', ...
                        'plot_correlograms.py failed:\n%s', result);
                end
                fprintf('%s', result);
                nOk = nOk + 1;
            catch ME
                nErr = nErr + 1;
                fprintf('       ERROR: %s\n', ME.message);
                errors{end+1} = sprintf('%s: %s', label, ME.message); %#ok<AGROW>
            end
        end
    end

    % --- Study-independent SI panels --------------------------------------
    if doSI
        calls = {
            'fig_mea_channel_map', {}
        };
        [nTotal, nOk, nErr, errors] = run_calls( ...
            calls, nTotal, nOk, nErr, errors);
    end

    % --- Summary ----------------------------------------------------------
    fprintf('\n==== run_figures summary ====\n');
    fprintf('  OK  : %d\n', nOk);
    fprintf('  FAIL: %d\n', nErr);
    fprintf('  Time: %.1f s\n', toc(t0));
    if nErr > 0
        fprintf('\nFailures:\n');
        for e = 1:numel(errors)
            fprintf('  - %s\n', errors{e});
        end
        error('run_figures:SomeFailed', '%d of %d calls failed.', nErr, nTotal);
    end
end


% =========================================================================
function [nTotal, nOk, nErr, errors] = run_calls(calls, nTotal, nOk, nErr, errors)
%RUN_CALLS Execute a cell array of {funcName, {args}} with error handling.
    for k = 1:size(calls, 1)
        fname = calls{k, 1};
        args  = calls{k, 2};
        nTotal = nTotal + 1;
        label  = sprintf('%s(%s)', fname, argstr(args));
        try
            fprintf('  [%2d] %s\n', nTotal, label);
            fcn = str2func(fname);
            fcn(args{:});
            nOk = nOk + 1;
        catch ME
            nErr = nErr + 1;
            fprintf('       ERROR: %s\n', ME.message);
            errors{end+1} = sprintf('%s: %s', label, ME.message); %#ok<AGROW>
            if ~isempty(ME.stack)
                fprintf('       at %s:%d\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
        close all;
    end
end

function s = argstr(args)
    parts = cell(1, numel(args));
    for i = 1:numel(args)
        if ischar(args{i}) || isstring(args{i})
            parts{i} = sprintf('''%s''', char(args{i}));
        elseif isstruct(args{i})
            parts{i} = '<struct>';
        elseif isnumeric(args{i}) || islogical(args{i})
            parts{i} = mat2str(args{i});
        else
            parts{i} = class(args{i});
        end
    end
    s = strjoin(parts, ',');
end
