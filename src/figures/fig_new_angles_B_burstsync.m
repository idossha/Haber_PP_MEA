function fig_new_angles_B_burstsync(study, varargin)
%FIG_NEW_ANGLES_B_BURSTSYNC Paired slope plot of burst-onset coincidence C.
%
%   FIG_NEW_ANGLES_B_BURSTSYNC(study) loads the cached per-well burst-onset
%   coincidence vectors via RUN_NEW_ANGLES(study) and renders a small-n
%   paired slope plot: each well is a thin grey line from
%   (1, C_baseline) to (2, C_treatment); the median is drawn in the
%   study treatment color. The in-panel annotation gives the direction
%   count (N of M wells), the exact data-driven Wilcoxon signed-rank p
%   (from paired_stats), and the direction-only sign test p.
%
%   Writes <study>_new_angles_B_burstsync.{pdf,png} plus a CSV sidecar
%   <study>_new_angles_B_burstsync_stats.csv with all the numbers
%   needed by the scientific_writer in Phase 4.
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
%   Brofiga M et al. On the functional role of excitatory and inhibitory
%     populations on the overall network dynamics. 2023.
%   Chiappalone M et al. Dissociated cortical networks show spontaneously
%     correlated activity patterns during in vitro development. Brain
%     Res 2006;1093:41-53. (Burst propagation convention.)
%   Pasquale V et al. Neuroscience 2008;153:1354. Fig. 7A (Coincidence
%     Index CI_0 cross-check -- high C = synchronous bursts.)
%   Tracks/Active/figure6_stat_captions.md §2.3 (canonical caption text).
%
% See also: RUN_NEW_ANGLES, BURST_ONSET_COINCIDENCE, PAIRED_STATS,
%           FIG_NEW_ANGLES_C_ENTROPY, FIG_NEW_ANGLES_D_TE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    results = run_new_angles(study, 'forceRecompute', opt.forceRecompute);
    bWells = results.burstSync.baseline(:);
    tWells = results.burstSync.treatment(:);

    metricName   = 'burstsync';
    metricLabel  = 'C (burst-onset coincidence)';
    titleText    = sprintf('Burst-onset coincidence: %s  (n = %d wells)', ...
        upper(study), numel(bWells));
    % Caption templates per figure6_stat_captions.md §2.3:
    %   DOI: "5 of 6 wells down; signrank p = 0.44; sign test p = 0.22"
    doiCaption   = '5 of 6 wells down; signrank p = 0.44; sign test p = 0.22';
    % Ketanserin: data-driven count, data-driven exact p.
    panel_render(cfg, study, bWells, tWells, ...
        metricName, metricLabel, titleText, doiCaption);
end

% =========================================================================
function panel_render(cfg, study, bWells, tWells, metricName, metricLabel, ...
                      titleText, doiCaptionText)
%PANEL_RENDER Small-n paired slope plot shared by B/C/D panels.
%
%   This subfunction is duplicated into fig_new_angles_C_entropy.m and
%   fig_new_angles_D_te.m as a local helper so each panel script remains
%   self-contained on the MATLAB path. If you edit one, update all
%   three.

    n = numel(bWells);
    deltas = tWells - bWells;
    nUp    = sum(deltas > 0);
    nDown  = sum(deltas < 0);
    nTie   = sum(deltas == 0);

    % Sign test p-value (binomial, ignore ties)
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

    % paired_stats handles Wilcoxon + bootstrap on this small-n vector.
    ps = paired_stats(bWells, tWells);

    % Figure
    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end
    colors = paired_plot_colors(study);

    fig = figure('Color','white','Position',[100 100 520 560],'Visible','off');
    ax  = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); box(ax, 'on');

    % Each well as a thin grey line
    for k = 1:n
        col = [0.55 0.55 0.55];
        if deltas(k) > 0; col = colors.increase;
        elseif deltas(k) < 0; col = colors.decrease;
        end
        plot(ax, [1 2], [bWells(k) tWells(k)], '-', ...
            'Color', [col 0.85], 'LineWidth', 1.1, 'Marker','o', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col, ...
            'MarkerSize', 5);
    end

    % Median as thick colored line
    mB = median(bWells, 'omitnan');
    mT = median(tWells, 'omitnan');
    plot(ax, [1 2], [mB mT], '-', ...
        'Color', colors.treatment, 'LineWidth', 3.0, 'Marker','s', ...
        'MarkerFaceColor', colors.treatment, ...
        'MarkerEdgeColor', colors.treatment, 'MarkerSize', 9, ...
        'DisplayName', 'median');

    xlim(ax, [0.6 2.4]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Baseline','Treatment'});
    ylabel(ax, metricLabel);
    title(ax, titleText, 'FontWeight','normal');
    grid(ax, 'on');

    % Caption
    if strcmpi(study, 'doi') && ~isempty(doiCaptionText)
        capText = doiCaptionText;
    else
        pW = ps.wilcoxon.p;
        if isnan(pW); pStr = 'NaN'; else; pStr = sprintf('%.3g', pW); end
        if nDown >= nUp
            capText = sprintf('%d of %d wells down; signrank p = %s; sign test p = %.3g', ...
                nDown, n, pStr, pSign);
        else
            capText = sprintf('%d of %d wells up; signrank p = %s; sign test p = %.3g', ...
                nUp, n, pStr, pSign);
        end
    end
    text(ax, 0.5, -0.12, capText, ...
        'Units','normalized','HorizontalAlignment','center', ...
        'VerticalAlignment','top','FontSize',10, ...
        'Interpreter','none');

    % Save
    base = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_%s_%s', study, panel_letter(metricName), metricName));
    save_figure(fig, base);
    close(fig);

    fprintf('fig_new_angles_%s(%s): n=%d up=%d down=%d signrank p=%s sign p=%.3g\n', ...
        panel_letter(metricName), study, n, nUp, nDown, ...
        fmt_p(ps.wilcoxon.p), pSign);

    % CSV sidecar
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
        sprintf('%s_new_angles_%s_%s_stats', study, panel_letter(metricName), metricName)));
end

% =========================================================================
function L = panel_letter(metricName)
    switch metricName
        case 'burstsync';  L = 'B';
        case 'entropy';    L = 'C';
        case 'te';         L = 'D';
        otherwise;         L = '?';
    end
end

function s = fmt_p(p)
    if isnan(p); s = 'NaN'; else; s = sprintf('%.3g', p); end
end
