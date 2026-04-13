function fig_rate_scatter(metric, study, varargin)
%FIG_RATE_SCATTER  Baseline vs treatment scatter plot for a rate metric.
%
%   FIG_RATE_SCATTER('spike', 'doi')  plots per-channel baseline firing rate
%   on the x-axis vs DOI firing rate on the y-axis (log-log). Points above
%   the identity line increased; points below decreased.
%
%   FIG_RATE_SCATTER('burst', 'doi')  same for burst rate.
%
% OUTPUTS:
%   PNG + PDF written to cfg.paths.figures_out.

    cfg = project_config();

    p = inputParser;
    addRequired(p, 'metric', @(s) any(strcmpi(s, {'spike','burst'})));
    addRequired(p, 'study',  @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels', cfg.channels.default);
    addParameter(p, 'ignoreSilentChannels', true);
    addParameter(p, 'silentMode', cfg.silent.mode);
    addParameter(p, 'minRateThreshold', []);
    parse(p, metric, study, varargin{:});
    opt = p.Results;
    metric = lower(opt.metric);
    study  = lower(opt.study);

    % Metric-specific config
    metricCfg = struct( ...
        'spike', struct('field','spikeRates', 'thresh',cfg.silent.min_rate_spike, ...
                        'label','Firing rate', 'unit','spikes min^{-1}'), ...
        'burst', struct('field','burstRates', 'thresh',cfg.silent.min_rate_burst, ...
                        'label','Burst rate',  'unit','bursts min^{-1}'));
    mc = metricCfg.(metric);
    if ~isempty(opt.minRateThreshold)
        mc.thresh = opt.minRateThreshold;
    end

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    [rB, rT] = load_pair_metric(pairs, opt.channels, mc.field, cfg);

    % Silent-channel exclusion
    if opt.ignoreSilentChannels
        active = rB >= mc.thresh | rT >= mc.thresh;
        rB = rB(active);
        rT = rT(active);
    end

    % Remove NaN
    valid = isfinite(rB) & isfinite(rT);
    rB = rB(valid);
    rT = rT(valid);
    nCh = numel(rB);

    increased = rT > rB;
    nInc = sum(increased);
    nDec = nCh - nInc;

    colors = paired_plot_colors(study);

    % --- Figure ---
    fig = figure('Visible', 'off', 'Position', [100 100 700 700], 'Color', 'w');
    set(fig, 'Units', 'inches', 'PaperUnits', 'inches', ...
        'PaperSize', [7 7], 'PaperPosition', [0 0 7 7]);
    hold on;

    % Floor for log scale: replace zeros with a small positive value.
    % When thresh is 0 (e.g. burst rate), fall back to a data-driven floor
    % so the log axis is well-defined.
    if mc.thresh > 0
        floor_val = mc.thresh * 0.5;
    else
        posVals = [rB(rB > 0); rT(rT > 0)];
        if ~isempty(posVals)
            floor_val = min(posVals) * 0.5;
        else
            floor_val = 0.1;
        end
    end
    rBplot = max(rB, floor_val);
    rTplot = max(rT, floor_val);

    % Axis range
    allVals = [rBplot; rTplot];
    axMin = min(allVals) * 0.7;
    axMax = max(allVals) * 1.4;

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([axMin axMax]);
    ylim([axMin axMax]);
    axis square;

    % Scatter: color each point by increase/decrease
    dotSize  = 50;
    dotAlpha = 0.8;
    pointColor = repmat(colors.decrease, nCh, 1);
    pointColor(increased, :) = repmat(colors.increase, nInc, 1);

    scatter(rBplot, rTplot, dotSize, pointColor, 'filled', ...
        'MarkerFaceAlpha', dotAlpha);

    % Dummy handles for legend (before identity line so line stays on top)
    hInc = scatter(nan, nan, dotSize, colors.increase, 'filled');
    hDec = scatter(nan, nan, dotSize, colors.decrease, 'filled');
    legend([hInc hDec], ...
        {sprintf('Increased (%d, %.0f%%)', nInc, 100*nInc/nCh), ...
         sprintf('Decreased (%d, %.0f%%)', nDec, 100*nDec/nCh)}, ...
        'Location', 'northwest', 'FontSize', 12, 'Box', 'off');

    xlabel(sprintf('Baseline %s (%s)', mc.label, mc.unit), ...
        'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    ylabel(sprintf('%s %s (%s)', labels.treatment, mc.label, mc.unit), ...
        'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 16, 'FontName', 'Arial', 'LineWidth', 1.5, ...
        'TickDir', 'out', 'Box', 'off');

    % Identity line (y = x) — drawn LAST so it sits on top of all scatter objects
    limRange = xlim;
    hLine = plot(limRange, limRange, 'k--', 'LineWidth', 1.5);
    uistack(hLine, 'top');
    set(hLine, 'HandleVisibility', 'off');

    % --- Save ---
    outBase = fullfile(cfg.paths.figures_out, ...
        sprintf('%s_%s_rate_scatter', study, metric));
    save_figure(fig, [outBase '.png']);
    exportgraphics(fig, [outBase '.pdf'], 'ContentType', 'vector');
    close(fig);

    fprintf('fig_rate_scatter(%s, %s): n=%d, increased=%d (%.1f%%), decreased=%d (%.1f%%)\n', ...
        metric, study, nCh, nInc, 100*nInc/nCh, nDec, 100*nDec/nCh);
    fprintf('  saved: %s.png\n', outBase);
end
