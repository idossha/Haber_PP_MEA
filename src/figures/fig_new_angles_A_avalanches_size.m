function fig_new_angles_A_avalanches_size(study, varargin)
%FIG_NEW_ANGLES_A_AVALANCHES_SIZE Avalanche size distribution log-log overlay.
%
%   FIG_NEW_ANGLES_A_AVALANCHES_SIZE(study) loads the cached per-pair
%   avalanche size histograms via RUN_NEW_ANGLES(study) and plots the
%   baseline-vs-treatment P(size) distributions on a single log-log
%   axis for the ONE representative pair that maximises
%   |alphaLS_treat - alphaLS_base|. A dashed reference line at slope
%   -1.5 (Pasquale 2008 p.1358; Massobrio 2015 p.2) is overlaid through
%   the baseline median for visual anchoring.
%
%   Writes <study>_new_angles_A_avalanches_size.{pdf,png} plus a
%   CSV sidecar <study>_new_angles_A_avalanches_size_stats.csv. The
%   CSV contains the representative-pair index and label and the LS /
%   MLE fit results for both baseline and treatment.
%
% INPUTS:
%   study  -  'doi' or 'ket' (case-insensitive).
%
% Name-value options:
%   'forceRecompute'  - forwarded to run_new_angles (default false).
%
% OUTPUTS:
%   PDF + PNG written to cfg.paths.figures_out via SAVE_FIGURE.
%   CSV sidecar via EXPORT_FIGURE_STATS.
%
% References:
%   Pasquale V, Massobrio P, Bologna LL, Chiappalone M, Martinoia S.
%     Self-organization and neuronal avalanches in networks of dissociated
%     cortical neurons. Neuroscience 2008;153:1354. p.1358, p.1361
%     (critical size exponent alpha = -1.5).
%   Massobrio P, Pasquale V, Martinoia S. Self-organized criticality in
%     cortical assemblies. Sci Rep 2015;5:10578.
%
% See also: RUN_NEW_ANGLES, FIG_NEW_ANGLES_E_AVALANCHE_CLASS.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'forceRecompute', false, @(x) islogical(x) && isscalar(x));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    %%% ----- load cached metrics -----------------------------------------
    results = run_new_angles(study, 'forceRecompute', opt.forceRecompute);
    nPairs = results.nPairs;
    if nPairs < 1
        error('fig_new_angles_A_avalanches_size:NoData', ...
            'No pairs available for study %s.', study);
    end

    %%% ----- pick representative pair by |delta alphaLS| -----------------
    deltas = nan(nPairs, 1);
    for k = 1:nPairs
        aB = results.avalanches.baseline(k).alphaLS;
        aT = results.avalanches.treatment(k).alphaLS;
        if isfinite(aB) && isfinite(aT)
            deltas(k) = abs(aT - aB);
        end
    end
    [maxDelta, repIdx] = max(deltas);
    if ~isfinite(maxDelta)
        % Fall back to pair 1 if no pair has finite fits.
        repIdx = 1;
    end
    repLabel = results.pairLabels{repIdx};

    bRow = results.avalanches.baseline(repIdx);
    tRow = results.avalanches.treatment(repIdx);

    %%% ----- extract size histograms -------------------------------------
    % sizeCounts is a 1 x nCh histogram (bin i = #avalanches of size i).
    bCounts = double(bRow.sizeCounts(:)');
    tCounts = double(tRow.sizeCounts(:)');
    bx = 1:numel(bCounts);
    tx = 1:numel(tCounts);

    % Normalise to probability for visual overlay
    bProb = bCounts / max(sum(bCounts), 1);
    tProb = tCounts / max(sum(tCounts), 1);

    % Drop zero bins so loglog is well defined
    bkeep = bProb > 0 & bx >= 1;
    tkeep = tProb > 0 & tx >= 1;

    %%% ----- figure ------------------------------------------------------
    if ~exist(cfg.paths.figures_out, 'dir')
        mkdir(cfg.paths.figures_out);
    end

    colors = paired_plot_colors(study);

    fig = figure('Color','white','Position',[100 100 620 520],'Visible','off');
    ax = axes(fig); %#ok<LAXES>
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    set(ax, 'XScale','log', 'YScale','log');

    plot(ax, bx(bkeep), bProb(bkeep), 'o-', ...
        'Color', colors.baseline,  'MarkerFaceColor', colors.baseline,  ...
        'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Baseline');
    plot(ax, tx(tkeep), tProb(tkeep), 's-', ...
        'Color', colors.treatment, 'MarkerFaceColor', colors.treatment, ...
        'LineWidth', 1.4, 'MarkerSize', 5, 'DisplayName', 'Treatment');

    % Slope = -1.5 reference line (Pasquale 2008 p.1358 critical exponent).
    % Anchor at the baseline median so the reference sits visually on the
    % observed curve.
    refSlope = -1.5;
    if any(bkeep)
        mid = find(bkeep, 1, 'first') + floor(sum(bkeep) / 2) - 1;
        if mid < 1 || mid > numel(bx); mid = find(bkeep, 1, 'first'); end
        x0 = bx(mid);  y0 = bProb(mid);
        xr = logspace(0, log10(max(bx(bkeep))), 40);
        yr = y0 * (xr / x0) .^ refSlope;
        plot(ax, xr, yr, '--', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.2, 'DisplayName', 'slope = -1.5 (Pasquale 2008)');
    end

    xlabel(ax, 'Avalanche size (# electrodes)');
    ylabel(ax, 'Probability P(size)');
    title(ax, sprintf('%s  |  %s\n\\alpha_{LS} base = %.2f, treat = %.2f', ...
        upper(study), repLabel, bRow.alphaLS, tRow.alphaLS), ...
        'Interpreter','tex','FontWeight','normal');
    legend(ax, 'Location','southwest','Box','off');

    % Annotate MLE values + classification in a text box.
    annot = sprintf(['\\alpha_{MLE} base = %.2f (KSp=%.2f) [%s]\n', ...
                     '\\alpha_{MLE} treat = %.2f (KSp=%.2f) [%s]'], ...
        bRow.alphaMLE, bRow.alphaMLE_ksp, bRow.classification, ...
        tRow.alphaMLE, tRow.alphaMLE_ksp, tRow.classification);
    text(ax, 0.98, 0.98, annot, ...
        'Units','normalized','HorizontalAlignment','right', ...
        'VerticalAlignment','top','Interpreter','tex','FontSize',9, ...
        'BackgroundColor',[1 1 1 0.75],'EdgeColor',[0.7 0.7 0.7]);

    %%% ----- save figure -------------------------------------------------
    base = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_A_avalanches_size', study));
    save_figure(fig, base);
    close(fig);

    fprintf('fig_new_angles_A_avalanches_size(%s): rep=pair %d |dA|=%.2f saved %s.{pdf,png}\n', ...
        study, repIdx, ifelse(isfinite(maxDelta), maxDelta, NaN), base);

    %%% ----- CSV sidecar -------------------------------------------------
    stats = struct( ...
        'study',                       study, ...
        'rep_pair_index',              repIdx, ...
        'rep_pair_label',              repLabel, ...
        'alpha_ls_baseline',           bRow.alphaLS, ...
        'alpha_ls_treatment',          tRow.alphaLS, ...
        'alpha_mle_baseline',          bRow.alphaMLE, ...
        'alpha_mle_treatment',         tRow.alphaMLE, ...
        'ks_p_baseline',               bRow.alphaMLE_ksp, ...
        'ks_p_treatment',              tRow.alphaMLE_ksp, ...
        'classification_baseline',     char(bRow.classification), ...
        'classification_treatment',    char(tRow.classification), ...
        'n_avalanches_baseline',       bRow.nAvalanches, ...
        'n_avalanches_treatment',      tRow.nAvalanches, ...
        'bin_ms_baseline',             bRow.binMs, ...
        'bin_ms_treatment',            tRow.binMs);
    export_figure_stats(stats, fullfile(cfg.paths.figures_out, ...
        sprintf('%s_new_angles_A_avalanches_size_stats', study)));
end

% =========================================================================
function v = ifelse(cond, a, b)
    if cond; v = a; else, v = b; end
end
