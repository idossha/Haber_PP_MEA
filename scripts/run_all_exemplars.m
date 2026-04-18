function run_all_exemplars(varargin)
%RUN_ALL_EXEMPLARS Generate traces, raster, and connectivity for every pair.
%
%   run_all_exemplars()                   DOI then ketanserin, all panels.
%   run_all_exemplars('study','doi')      DOI only.
%   run_all_exemplars('connectivity', false)  skip connectivity panels.
%
%   Output directory: <repo>/output/fig{2,4}/exemplars/
%   Naming:  {study}_pair{k}_traces.png/pdf
%            {study}_pair{k}_raster.png/pdf
%            {study}_pair{k}_connectivity.png/pdf + _stats.json/csv
%
%   After inspection, select which pair to feature in main Figures 2-5.
%   All per-pair panels are included in the Supplementary Information
%   (Figures S3-S8).
%
% See also: FIG_RATE_PANEL_TRACES, FIG_RATE_PANEL_RASTER,
%           FIG_CONNECTIVITY_EXEMPLAR, RUN_FIGURES.

    p = inputParser;
    addParameter(p, 'study',        'all', @(s) any(strcmpi(s, {'all','doi','ket'})));
    addParameter(p, 'connectivity', true,  @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    cfg = project_config();

    switch lower(opt.study)
        case 'all'; studies = {'doi','ket'};
        case 'doi'; studies = {'doi'};
        case 'ket'; studies = {'ket'};
    end

    t0 = tic;
    nOk = 0; nErr = 0; errors = {};

    fprintf('\n==== run_all_exemplars ====\n');
    fprintf('studies     : %s\n', strjoin(studies, ', '));
    fprintf('connectivity: %d\n\n', opt.connectivity);

    for s = 1:numel(studies)
        study = studies{s};
        [pairs, ~] = get_pairs_and_labels(cfg, study);
        nPairs = numel(pairs);

        rateExDir = output_path(cfg, study, 'rates', 'exemplars');
        connExDir = output_path(cfg, study, 'connectivity', 'exemplars');
        connStDir = output_path(cfg, study, 'connectivity', 'stats');
        if ~exist(rateExDir, 'dir'); mkdir(rateExDir); end

        fprintf('--- %s: %d pairs ---\n', upper(study), nPairs);
        fprintf('    rate exemplars : %s\n', rateExDir);

        for k = 1:nPairs
            tag = sprintf('%s_pair%d', study, k);
            fprintf('\n  [pair %d/%d] %s\n', k, nPairs, tag);
            fprintf('    baseline : %s\n', pairs(k).baseline);
            fprintf('    treatment: %s\n', pairs(k).treatment);

            % --- Traces ---
            try
                fprintf('    traces...');
                fig_rate_panel_traces(study, ...
                    'pairIndex', k, ...
                    'outName',   [tag '_traces'], ...
                    'outDir',    rateExDir);
                nOk = nOk + 1;
                fprintf(' OK\n');
            catch ME
                nErr = nErr + 1;
                errors{end+1} = sprintf('%s traces: %s', tag, ME.message); %#ok<AGROW>
                fprintf(' ERROR: %s\n', ME.message);
            end
            close all;

            % --- Raster ---
            try
                fprintf('    raster...');
                fig_rate_panel_raster(study, ...
                    'pairIndex', k, ...
                    'outName',   [tag '_raster'], ...
                    'outDir',    rateExDir);
                nOk = nOk + 1;
                fprintf(' OK\n');
            catch ME
                nErr = nErr + 1;
                errors{end+1} = sprintf('%s raster: %s', tag, ME.message); %#ok<AGROW>
                fprintf(' ERROR: %s\n', ME.message);
            end
            close all;

            % --- Connectivity ---
            if opt.connectivity
                if ~exist(connExDir, 'dir'); mkdir(connExDir); end
                if ~exist(connStDir, 'dir'); mkdir(connStDir); end
                try
                    fprintf('    connectivity...');
                    fig_connectivity_exemplar(study, ...
                        'pairIndex', k, ...
                        'outName',   [tag '_connectivity'], ...
                        'outDir',    connExDir, ...
                        'statsDir',  connStDir);
                    nOk = nOk + 1;
                    fprintf(' OK\n');
                catch ME
                    nErr = nErr + 1;
                    errors{end+1} = sprintf('%s connectivity: %s', tag, ME.message); %#ok<AGROW>
                    fprintf(' ERROR: %s\n', ME.message);
                end
                close all;
            end
        end
    end

    fprintf('\n==== run_all_exemplars summary ====\n');
    fprintf('  OK   : %d\n', nOk);
    fprintf('  FAIL : %d\n', nErr);
    fprintf('  Time : %.1f s\n', toc(t0));
    if nErr > 0
        fprintf('\nFailures:\n');
        for e = 1:numel(errors)
            fprintf('  - %s\n', errors{e});
        end
        error('run_all_exemplars:SomeFailed', '%d of %d panels failed.', nErr, nOk + nErr);
    end
end
