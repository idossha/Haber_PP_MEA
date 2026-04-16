function run_figures(varargin)
%RUN_FIGURES Batch entrypoint that renders every figure for a study.
%
%   run_figures()                  renders DOI then ketanserin.
%   run_figures('study','doi')     renders DOI only.
%   run_figures('study','ket')     renders ketanserin only.
%   run_figures('connectivity',    false)  skips connectivity panels.
%
%   Preconditions: pipeline/preprocess_and_save.m has been run so the
%   v2.0 caches exist under <repo>/cache/. If a cache is missing the
%   figure script will raise :MissingCache and this driver will log the
%   error but continue with the remaining figures.
%
%   Output: PNG files under <repo>/figures_out/ (see
%   cfg.paths.figures_out in src/config/project_config.m).

    p = inputParser;
    addParameter(p, 'study',        'all', @(s) any(strcmpi(s, {'all','doi','ket'})));
    addParameter(p, 'connectivity', true,  @(x) islogical(x) && isscalar(x));
    addParameter(p, 'newAngles',    true,  @(x) islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    studies = {};
    switch lower(opt.study)
        case 'all'; studies = {'doi','ket'};
        case 'doi'; studies = {'doi'};
        case 'ket'; studies = {'ket'};
    end

    t0 = tic;
    fprintf('\n==== run_figures ====\n');
    fprintf('studies     : %s\n', strjoin(studies, ', '));
    fprintf('connectivity: %d\n', opt.connectivity);
    fprintf('newAngles   : %d\n\n', opt.newAngles);

    nTotal = 0; nOk = 0; nErr = 0;
    errors = {};

    for s = 1:numel(studies)
        study = studies{s};
        calls = {
            'fig_rate_panel_traces',          {study};
            'fig_rate_panel_raster',          {study};
            'fig_rate_scatter',               {'spike', study};
            'fig_rate_scatter',               {'burst', study};
            'fig_spike_rate',                 {study};
            'fig_burst_rate',                 {study};
            'fig_spike_pct_change_violin',    {study};
            'fig_burst_pct_change_violin',    {study};
            'fig_rate_change_bar',            {'spike', study};
            'fig_rate_change_bar',            {'burst', study};
            'fig_stats_bootstrap',            {study};
        };
        if opt.connectivity
            calls = [calls; { ...
                'fig_connectivity_summary',   {study};
                'fig_connectivity_exemplar',  {study}}];
        end
        if opt.newAngles
            calls = [calls; { ...
                'fig_new_angles_A_avalanches_size', {study};
                'fig_new_angles_B_burstsync',       {study};
                'fig_new_angles_C_entropy',         {study};
                'fig_new_angles_D_te',              {study};
                'fig_new_angles_E_avalanche_class', {study}}];
        end

        for k = 1:size(calls, 1)
            fname = calls{k, 1};
            args  = calls{k, 2};
            nTotal = nTotal + 1;
            label = sprintf('%s(%s)', fname, argstr(args));
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

    fprintf('\n==== run_figures summary ====\n');
    fprintf('  OK  : %d\n', nOk);
    fprintf('  FAIL: %d\n', nErr);
    fprintf('  Total time: %.1f s\n', toc(t0));
    if nErr > 0
        fprintf('\nFailures:\n');
        for e = 1:numel(errors)
            fprintf('  - %s\n', errors{e});
        end
        error('run_figures:SomeFailed', '%d of %d figures failed.', nErr, nTotal);
    end
end

% -------------------------------------------------------------------------
function s = argstr(args)
    parts = cell(1, numel(args));
    for i = 1:numel(args)
        if ischar(args{i}) || isstring(args{i})
            parts{i} = sprintf('''%s''', char(args{i}));
        else
            parts{i} = mat2str(args{i});
        end
    end
    s = strjoin(parts, ',');
end
