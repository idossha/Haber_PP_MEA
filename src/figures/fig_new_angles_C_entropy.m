function fig_new_angles_C_entropy(study, varargin)
%FIG_NEW_ANGLES_C_ENTROPY Paired slope plot of per-well Shannon entropy H.
%
%   FIG_NEW_ANGLES_C_ENTROPY(study) loads the cached per-well median
%   Shannon entropy H (50 ms bin, clip-3 alphabet, Varley 2024 §2.2)
%   via RUN_NEW_ANGLES(study) and renders a small-n paired slope plot:
%   each well is a thin grey line from (1, H_baseline) to
%   (2, H_treatment); the median is drawn in the study treatment color.
%   The in-panel annotation gives the direction count (N of M wells),
%   the exact data-driven Wilcoxon signed-rank p, and the
%   direction-only sign test p.
%
%   Writes <study>_new_angles_C_entropy.{pdf,png} plus a CSV sidecar
%   <study>_new_angles_C_entropy_stats.csv.
%
% INPUTS:
%   study  -  'doi' or 'ket' (case-insensitive).
%
% Name-value options:
%   'forceRecompute'  - forwarded to run_new_angles (default false).
%
% OUTPUTS:
%   PDF + PNG via SAVE_FIGURE; CSV sidecar via EXPORT_FIGURE_STATS.
%
% References:
%   Varley TF, Pope M, Grazzi LF, Pittman-Polletta B, Sporns O.
%     Information processing dynamics in neural networks of macroscopic
%     brain organoids. J Neural Eng 2024. (§2.2 -- 50 ms bin, clip-3
%     Shannon H convention; p.1425 log-H increase under DPT.)
%   Carhart-Harris RL et al. The entropic brain: a theory of conscious
%     states informed by neuroimaging research with psychedelic drugs.
%     Front Hum Neurosci 2014.
%   Tracks/Active/figure6_stat_captions.md §2.2 (canonical caption text).
%
% See also: RUN_NEW_ANGLES, SPIKE_ENTROPY_SHANNON, PAIRED_STATS,
%           FIG_NEW_ANGLES_B_BURSTSYNC, FIG_NEW_ANGLES_D_TE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    results = run_new_angles(study, 'forceRecompute', opt.forceRecompute);
    bWells = results.entropy.baseline(:);
    tWells = results.entropy.treatment(:);

    metricName   = 'entropy';
    metricLabel  = 'Shannon entropy H (bits)';
    titleText    = sprintf('Spike-train Shannon entropy: %s  (n = %d wells)', ...
        upper(study), numel(bWells));
    % Caption per figure6_stat_captions.md §2.2 (DOI):
    %   "4 of 6 wells up; signrank p = 0.44; sign test p = 0.69"
    doiCaption   = '4 of 6 wells up; signrank p = 0.44; sign test p = 0.69';

    panel_render(cfg, study, bWells, tWells, ...
        metricName, metricLabel, titleText, doiCaption);
end

% =========================================================================
function panel_render(cfg, study, bWells, tWells, metricName, metricLabel, ...
                      titleText, doiCaptionText)
%PANEL_RENDER Small-n paired slope plot (duplicated across B/C/D panels).

    n = numel(bWells);
    deltas = tWells - bWells;
    nUp    = sum(deltas > 0);
    nDown  = sum(deltas < 0);
    nTie   = sum(deltas == 0);

    nDir = nUp + nDown;
    if nDir > 0
        k = min(nUp, nDown);
        pSign = 0;
        for j = 0:k
            pSign = pSign + nchoosek(nDir, j) * 0.5^nDir;
        end
        pSign = min(2 * pSign, 1);
    else
        pSign = NaN;
    end

    ps = paired_stats(bWells, tWells);

    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end
    colors = paired_plot_colors(study);

    fig = figure('Color','white','Position',[100 100 520 560],'Visible','off');
    ax  = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); box(ax, 'on');

    for k = 1:n
        col = [0.55 0.55 0.55];
        if deltas(k) > 0; col = colors.increase;
        elseif deltas(k) < 0; col = colors.decrease;
        end
        plot(ax, [1 2], [bWells(k) tWells(k)], '-', ...
            'Color', [col 0.85], 'LineWidth', 1.1, 'Marker','o', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col, 'MarkerSize', 5);
    end

    mB = median(bWells, 'omitnan');
    mT = median(tWells, 'omitnan');
    plot(ax, [1 2], [mB mT], '-', ...
        'Color', colors.treatment, 'LineWidth', 3.0, 'Marker','s', ...
        'MarkerFaceColor', colors.treatment, ...
        'MarkerEdgeColor', colors.treatment, 'MarkerSize', 9);

    xlim(ax, [0.6 2.4]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Baseline','Treatment'});
    ylabel(ax, metricLabel);
    title(ax, titleText, 'FontWeight','normal');
    grid(ax, 'on');

    if strcmpi(study, 'doi') && ~isempty(doiCaptionText)
        capText = doiCaptionText;
    else
        pW = ps.wilcoxon.p;
        if isnan(pW); pStr = 'NaN'; else; pStr = sprintf('%.3g', pW); end
        if nUp >= nDown
            capText = sprintf('%d of %d wells up; signrank p = %s; sign test p = %.3g', ...
                nUp, n, pStr, pSign);
        else
            capText = sprintf('%d of %d wells down; signrank p = %s; sign test p = %.3g', ...
                nDown, n, pStr, pSign);
        end
    end
    text(ax, 0.5, -0.12, capText, ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',10,'Interpreter','none');

    base = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_C_%s', study, metricName));
    save_figure(fig, base);
    close(fig);

    fprintf('fig_new_angles_C_entropy(%s): n=%d up=%d down=%d signrank p=%s sign p=%.3g\n', ...
        study, n, nUp, nDown, fmt_p(ps.wilcoxon.p), pSign);

    stats = struct( ...
        'study',             study, ...
        'metric',            metricName, ...
        'n_pairs',           n, ...
        'n_up',              nUp, ...
        'n_down',            nDown, ...
        'n_tie',             nTie, ...
        'median_baseline',   mB, ...
        'median_treatment',  mT, ...
        'median_delta',      ps.medianDelta, ...
        'ci_delta_lo',       ps.bootstrap.ciDelta(1), ...
        'ci_delta_hi',       ps.bootstrap.ciDelta(2), ...
        'p_bootstrap',       ps.bootstrap.pDelta, ...
        'p_wilcoxon',        ps.wilcoxon.p, ...
        'p_sign_test',       pSign, ...
        'hedges_g_av',       ps.hedgesGav, ...
        'caption_text',      capText);
    export_figure_stats(stats, fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_C_%s_stats', study, metricName)));
end

% =========================================================================
function s = fmt_p(p)
    if isnan(p); s = 'NaN'; else; s = sprintf('%.3g', p); end
end
