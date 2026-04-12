function fig_new_angles_D_te(study, varargin)
%FIG_NEW_ANGLES_D_TE Paired slope plot of mean delay-1 transfer entropy.
%
%   FIG_NEW_ANGLES_D_TE(study) loads the cached per-well top-decile
%   mean delay-1 transfer entropy (Ito 2011 p.4 D1TE) via
%   RUN_NEW_ANGLES(study) and renders a small-n paired slope plot:
%   each well is a thin grey line from (1, TE_baseline) to
%   (2, TE_treatment); the median is drawn in the study treatment color.
%
%   The caption prints the data-driven exact Wilcoxon signed-rank
%   p-value (for DOI: **p = 0.0625**, NOT 0.094, per
%   Tracks/Active/figure6_stat_captions.md §2.1) and flags the
%   direction disagreement with Varley 2024 via a footnote pointer to
%   Discussion §"Information flow and criticality".
%
%   Writes <study>_new_angles_D_te.{pdf,png} plus a CSV sidecar
%   <study>_new_angles_D_te_stats.csv.
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
%   Ito S, Hansen ME, Heiland R, Lumsdaine A, Litke AM, Beggs JM.
%     Extending transfer entropy improves identification of effective
%     connectivity in a spiking cortical network model. PLoS ONE
%     2011;6:e27431. (p.4 D1TE definition.)
%   Varley TF et al. Information processing dynamics in neural networks
%     of macroscopic brain organoids. J Neural Eng 2024. (§p.1425
%     mTE decrease under DPT -- OPPOSITE sign; reconciliation in
%     Discussion.)
%   Tracks/Active/figure6_stat_captions.md §2.1 (canonical caption text,
%     exact p = 0.0625).
%   Tracks/Active/figure6_biological_interpretation.md §1 (TE vs Varley
%     reconciliation paragraph).
%
% See also: RUN_NEW_ANGLES, TRANSFER_ENTROPY_D1, PAIRED_STATS,
%           FIG_NEW_ANGLES_B_BURSTSYNC, FIG_NEW_ANGLES_C_ENTROPY.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    results = run_new_angles(study, 'forceRecompute', opt.forceRecompute);
    bWells = results.te.baseline(:);
    tWells = results.te.treatment(:);

    metricName   = 'te';
    metricLabel  = 'Mean top-decile TE (bits)';
    titleText    = sprintf('Delay-1 transfer entropy: %s  (n = %d wells)', ...
        upper(study), numel(bWells));
    % Caption per figure6_stat_captions.md §2.1 (DOI):
    %   "5 of 6 wells up; signrank p = 0.0625 (exact); sign test p = 0.22"
    %   asterisk footnote -> Discussion "Information flow and criticality"
    %   for Varley 2024 sign-reconciliation.
    doiCaption   = ['5 of 6 wells up; signrank p = 0.0625 (exact); ', ...
                    'sign test p = 0.22  *see Discussion'];

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

    % Footnote pointer to Discussion (Varley 2024 reconciliation)
    text(ax, 0.5, -0.18, '*reconciliation with Varley 2024 mTE: Discussion §Information flow and criticality', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',8,'FontAngle','italic', ...
        'Color',[0.35 0.35 0.35],'Interpreter','none');

    base = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_D_%s', study, metricName));
    save_figure(fig, base);
    close(fig);

    fprintf('fig_new_angles_D_te(%s): n=%d up=%d down=%d signrank p=%s sign p=%.3g\n', ...
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
        sprintf('%s_new_angles_D_%s_stats', study, metricName)));
end

% =========================================================================
function s = fmt_p(p)
    if isnan(p); s = 'NaN'; else; s = sprintf('%.3g', p); end
end
