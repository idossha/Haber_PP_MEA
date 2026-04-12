function fig_new_angles_E_avalanche_class(study, varargin)
%FIG_NEW_ANGLES_E_AVALANCHE_CLASS Per-recording avalanche classification table.
%
%   FIG_NEW_ANGLES_E_AVALANCHE_CLASS(study) loads the cached per-pair
%   avalanche classifications via RUN_NEW_ANGLES(study) and renders a
%   compact 2 x nPairs colored table (row 1 = baseline, row 2 =
%   treatment) showing sub / critical / super / indeterminate for each
%   recording. Each cell is annotated with a single-letter code
%   (S = subcritical, C = critical, X = supercritical, ? =
%   indeterminate) inside a categorically-coloured patch.
%
%   Writes <study>_new_angles_E_avalanche_class.{pdf,png} plus a
%   per-recording CSV with the full Clauset/LRT/Sethna numbers.
%
%   Classification tree anchored in Tracks/Active/phase1_decisions.md
%   §1 and avalanche_method_spec.md §10; colour palette is a
%   categorical (non-monotone) scheme so readers do not accidentally
%   read a gradient.
%
% INPUTS:
%   study  -  'doi' or 'ket' (case-insensitive).
%
% Name-value options:
%   'forceRecompute'  - forwarded to run_new_angles (default false).
%
% OUTPUTS:
%   PDF + PNG via SAVE_FIGURE; per-cell CSV sidecar written row-by-row
%   as a multi-line table (NOT via export_figure_stats, which is
%   single-row-oriented).
%
% References:
%   Pasquale V et al. Self-organization and neuronal avalanches in
%     networks of dissociated cortical neurons. Neuroscience
%     2008;153:1354. (Sub/critical/super shape rules.)
%   Massobrio P, Pasquale V, Martinoia S. Self-organized criticality in
%     cortical assemblies. Sci Rep 2015;5:10578. (Clauset MLE + KS
%     threshold 0.10; LRT alternatives set.)
%
% See also: RUN_NEW_ANGLES, AVALANCHES, FIG_NEW_ANGLES_A_AVALANCHES_SIZE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    results = run_new_angles(study, 'forceRecompute', opt.forceRecompute);
    nPairs = results.nPairs;
    if nPairs < 1
        error('fig_new_angles_E_avalanche_class:NoData', ...
            'No pairs available for study %s.', study);
    end

    %%% ----- build 2 x nPairs classification grid -----------------------
    classGrid = cell(2, nPairs);
    for k = 1:nPairs
        classGrid{1, k} = char(results.avalanches.baseline(k).classification);
        classGrid{2, k} = char(results.avalanches.treatment(k).classification);
    end

    % Categorical colour palette (colourblind-safe where possible)
    pal.subcritical   = [0.1725 0.4824 0.7137];   % #2C7BB6
    pal.critical      = [0.9922 0.6824 0.3804];   % #FDAE61
    pal.supercritical = [0.8431 0.0980 0.1098];   % #D7191C
    pal.indeterminate = [0.7412 0.7412 0.7412];   % #BDBDBD

    % Numeric matrix for imagesc (just used as a discrete colour index)
    codeMat = zeros(2, nPairs);
    letters = cell(2, nPairs);
    faceColors = zeros(2, nPairs, 3);
    for r = 1:2
        for c = 1:nPairs
            [codeMat(r, c), letters{r, c}, faceColors(r, c, :)] = ...
                class_to_code(classGrid{r, c}, pal);
        end
    end

    %%% ----- figure -----------------------------------------------------
    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end

    fig = figure('Color','white', ...
        'Position',[100 100 max(520, 80*nPairs + 180) 280], ...
        'Visible','off');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on');

    % Draw each cell as a coloured rectangle with the letter centered.
    for r = 1:2
        for c = 1:nPairs
            rgb = squeeze(faceColors(r, c, :))';
            rectangle(ax, 'Position', [c-0.5, (2-r)+0.0, 1, 1], ...
                'FaceColor', rgb, 'EdgeColor', [1 1 1], 'LineWidth', 1.5);
            textColor = best_contrast(rgb);
            text(ax, c, (2-r)+0.5, letters{r, c}, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontSize', 14, 'FontWeight','bold', ...
                'Color', textColor);
        end
    end

    xlim(ax, [0.5, nPairs + 0.5]);
    ylim(ax, [-0.1, 2.1]);
    set(ax, 'XTick', 1:nPairs, ...
            'XTickLabel', arrayfun(@(k) sprintf('p%d', k), 1:nPairs, 'UniformOutput', false), ...
            'YTick', [0.5 1.5], ...
            'YTickLabel', {'Treatment','Baseline'});
    set(ax, 'TickLength', [0 0], 'XAxisLocation','bottom');
    title(ax, sprintf('%s  |  Avalanche classification (n = %d wells)', ...
        upper(study), nPairs), 'FontWeight','normal');

    % Legend as discrete colour patches below the plot.
    legendLabels = {'S subcritical','C critical','X supercritical','? indeterminate'};
    legendColors = [pal.subcritical; pal.critical; pal.supercritical; pal.indeterminate];
    for k = 1:4
        patchX = 0.5 + (k-1)*(nPairs/4);
        patchY = -0.9;
        rectangle(ax, 'Position', [patchX, patchY, 0.5, 0.35], ...
            'FaceColor', legendColors(k, :), 'EdgeColor', [0.3 0.3 0.3]);
        text(ax, patchX + 0.6, patchY + 0.17, legendLabels{k}, ...
            'FontSize', 9, 'VerticalAlignment','middle');
    end
    ylim(ax, [-1.1, 2.1]);

    %%% ----- save figure ------------------------------------------------
    base = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_E_avalanche_class', study));
    save_figure(fig, base);
    close(fig);

    fprintf('fig_new_angles_E_avalanche_class(%s): n=%d pairs saved %s.{pdf,png}\n', ...
        study, nPairs, base);

    %%% ----- CSV sidecar (multi-row) ------------------------------------
    % export_figure_stats is single-row. Write the per-recording table
    % directly here so the scientific_writer (Phase 4) can pull all
    % Clauset / LRT / Sethna numbers for the SI Table §S3.
    csvPath = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_E_avalanche_class_stats.csv', study));
    fid = fopen(csvPath, 'w');
    if fid == -1
        warning('fig_new_angles_E_avalanche_class:CsvFailed', ...
            'Could not open %s for writing.', csvPath);
        return;
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, ['study,pair_idx,condition,classification,alpha_ls,', ...
        'alpha_mle,ks_p,lrt_preferred_model,sethna_residual,', ...
        'tail_mass_frac,n_avalanches\n']);
    for k = 1:nPairs
        write_av_row(fid, study, k, 'baseline',  results.avalanches.baseline(k));
        write_av_row(fid, study, k, 'treatment', results.avalanches.treatment(k));
    end
    clear cleanup;
end

% =========================================================================
function [code, letter, rgb] = class_to_code(classLabel, pal)
    switch lower(strtrim(classLabel))
        case 'subcritical';    code = 1; letter = 'S'; rgb = pal.subcritical;
        case 'critical';       code = 2; letter = 'C'; rgb = pal.critical;
        case 'supercritical';  code = 3; letter = 'X'; rgb = pal.supercritical;
        otherwise;             code = 4; letter = '?'; rgb = pal.indeterminate;
    end
end

% =========================================================================
function col = best_contrast(rgb)
    % Choose black or white text for max contrast on the given face
    % colour using the standard relative-luminance formula.
    lum = 0.2126*rgb(1) + 0.7152*rgb(2) + 0.0722*rgb(3);
    if lum > 0.55
        col = [0 0 0];
    else
        col = [1 1 1];
    end
end

% =========================================================================
function write_av_row(fid, study, pairIdx, cond, row)
    fprintf(fid, '%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d\n', ...
        study, pairIdx, cond, ...
        safe_str(row.classification), ...
        num(row.alphaLS), ...
        num(row.alphaMLE), ...
        num(row.alphaMLE_ksp), ...
        safe_str(row.lrtPreferredModel), ...
        num(row.sethnaResidual), ...
        num(row.tailMassFrac), ...
        row.nAvalanches);
end

function s = safe_str(v)
    if ischar(v); s = v;
    elseif isstring(v); s = char(v);
    else; s = '';
    end
end

function s = num(v)
    if isempty(v) || ~isscalar(v) || isnan(v)
        s = 'NaN';
    elseif isinf(v)
        if v > 0; s = 'Inf'; else; s = '-Inf'; end
    else
        s = sprintf('%.6g', v);
    end
end
