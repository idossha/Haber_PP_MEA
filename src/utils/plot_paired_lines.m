function plot_paired_lines(yBaseline, yTreatment, colors, statsString, outFile, varargin)
%PLOT_PAIRED_LINES Paired-line + half-violin figure for baseline vs treatment.
%
%   PLOT_PAIRED_LINES(yBaseline, yTreatment, colors, statsString, outFile)
%   renders the canonical paired-line figure used by fig_spike_rate /
%   fig_burst_rate: each baseline-treatment pair is connected by a colored
%   line (green for increase, red for decrease), translucent dots are drawn
%   on each side, half-violin densities are drawn behind the dots, and a
%   black IQR marker plus mean dot is overlaid (using the Tukey 1.5*IQR
%   inlier rule).
%
%   PLOT_PAIRED_LINES(..., 'Name', value, ...) with any of:
%     'title'          figure title (default '')
%     'yLabel'         y-axis label (default '')
%     'xTickLabels'    2x1 cell with baseline/treatment labels
%                      (default {'Baseline', 'Treatment'})
%     'yLimitMode'     'auto' (default) | 'robust'
%                      'robust' caps the y-axis at max(data cap,
%                      97.5th percentile * 1.3) so a handful of
%                      outliers don't squash the boxes.
%
% INPUTS:
%   yBaseline    -  Column vector of baseline values (one per channel).
%   yTreatment   -  Matching column vector of treatment values.
%   colors       -  Struct from paired_plot_colors().
%   statsString  -  Annotation string drawn in the top-right corner.
%   outFile      -  Absolute path for the saved PNG.
%
% OUTPUTS:
%   PNG file written to outFile (300 dpi).

    p = inputParser;
    addRequired(p,  'yBaseline');
    addRequired(p,  'yTreatment');
    addRequired(p,  'colors');
    addRequired(p,  'statsString');
    addRequired(p,  'outFile');
    addParameter(p, 'title',        '');
    addParameter(p, 'yLabel',       '');
    addParameter(p, 'xTickLabels',  {'Baseline', 'Treatment'});
    addParameter(p, 'yLimitMode',   'robust', @(s) any(strcmpi(s, {'auto','robust'})));
    addParameter(p, 'yScale',       'linear', @(s) any(strcmpi(s, {'linear','log'})));
    addParameter(p, 'yFloor',       0.3,      @(x) isscalar(x) && x > 0);
    parse(p, yBaseline, yTreatment, colors, statsString, outFile, varargin{:});
    opt = p.Results;

    nCh = numel(yBaseline);
    xBaseline  = ones(nCh, 1);
    xTreatment = repmat(1.22, nCh, 1);

    % For log-scale Y, clamp non-positive values to the floor so they
    % plot at the bottom of the visible range instead of -Inf. Keep the
    % unclamped originals for legend counts / stats annotation.
    if strcmpi(opt.yScale, 'log')
        yBasePlot  = max(yBaseline,  opt.yFloor);
        yTreatPlot = max(yTreatment, opt.yFloor);
    else
        yBasePlot  = yBaseline;
        yTreatPlot = yTreatment;
    end

    fig = figure('Visible', 'off', 'Position', [100, 100, 720, 720], 'Color', 'w');
    set(fig, 'Units', 'inches', 'PaperUnits', 'inches', ...
        'PaperSize', [7.2 7.2], 'PaperPosition', [0 0 7.2 7.2]);
    hold on;

    % Half-violin densities behind dots. For log scale, compute the KDE
    % in log space on the clamped data so the shape reflects the log
    % distribution.
    yAll  = [yBaseline; yTreatment];
    if strcmpi(opt.yScale, 'log')
        yAllPlot = [yBasePlot; yTreatPlot];
        yGrid = logspace(log10(min(yAllPlot)), log10(max(yAllPlot) + eps), 120);
    else
        yGrid = linspace(max(0, min(yAll) - 1), max(yAll) + 1, 100);
    end
    [fB, ~] = ksdensity(yBasePlot,  yGrid);
    [fT, ~] = ksdensity(yTreatPlot, yGrid);
    violinW = 0.06;

    if max(fB) > 0
        fB = (fB / max(fB)) * violinW;
        yB = yGrid(:);
        hpB = patch([0.98 - fB(:); 0.98 * ones(numel(yB), 1)], ...
                    [yB; flipud(yB)], colors.baseline, ...
                    'FaceAlpha', 0.6, 'EdgeColor', 'none');
        uistack(hpB, 'bottom');
        plot(0.98 - fB, yB, 'Color', colors.baselineOutline, 'LineWidth', 1.2);
    end
    if max(fT) > 0
        fT = (fT / max(fT)) * violinW;
        yT = yGrid(:);
        hpT = patch([1.24 * ones(numel(yT), 1); 1.24 + flipud(fT(:))], ...
                    [yT; flipud(yT)], colors.treatment, ...
                    'FaceAlpha', 0.6, 'EdgeColor', 'none');
        uistack(hpT, 'bottom');
        plot(1.24 + fT, yT, 'Color', colors.treatmentOutline, 'LineWidth', 1.2);
    end

    xSummaryBaseline  = 0.98 - 0.006;
    xSummaryTreatment = 1.24 + 0.006;

    % Connecting lines (store handles to build a legend).
    hIncrease = gobjects(0);
    hDecrease = gobjects(0);
    for i = 1:nCh
        if yTreatment(i) > yBaseline(i)
            lc = colors.increase;
            hLn = plot([xBaseline(i), xTreatment(i)], ...
                       [yBasePlot(i), yTreatPlot(i)], ...
                       'Color', lc, 'LineWidth', 1.2);
            if isempty(hIncrease); hIncrease = hLn; end
        else
            lc = colors.decrease;
            hLn = plot([xBaseline(i), xTreatment(i)], ...
                       [yBasePlot(i), yTreatPlot(i)], ...
                       'Color', lc, 'LineWidth', 1.2);
            if isempty(hDecrease); hDecrease = hLn; end
        end
    end

    scatter(xBaseline,  yBasePlot,  76, colors.baseline,  'filled', ...
        'MarkerFaceAlpha', colors.alphaDots);
    scatter(xTreatment, yTreatPlot, 76, colors.treatment, 'filled', ...
        'MarkerFaceAlpha', colors.alphaDots);

    plot_mean_iqr_marker(gca, xSummaryBaseline,  yBasePlot,  0.0035, 0.85, 5);
    plot_mean_iqr_marker(gca, xSummaryTreatment, yTreatPlot, 0.0035, 0.85, 5);

    % --- Axis range selection --------------------------------------------
    xlim([0.88, 1.30]);
    xticks([1, 1.22]);
    xticklabels(opt.xTickLabels);

    ymax = robust_ymax(yAll, opt.yLimitMode);
    if strcmpi(opt.yScale, 'log')
        % Floor so that log scale is well defined; anything below yFloor
        % is drawn AT yFloor (these are 'silent' channels under a
        % user-specified threshold).
        set(gca, 'YScale', 'log');
        yTop = ymax * 1.20;
        if yTop <= opt.yFloor
            yTop = opt.yFloor * 10;
        end
        ylim([opt.yFloor, yTop]);
        % Decade ticks (10^k) between floor and top.
        k0 = floor(log10(opt.yFloor));
        k1 = ceil(log10(yTop));
        yticks(10 .^ (k0:k1));
    else
        yTop = ymax * 1.05;
        ylim([0, yTop]);
        yTickStep = nice_tick_step(yTop, 8);
        yticks(0:yTickStep:ceil(yTop / yTickStep) * yTickStep);
    end

    if ~isempty(opt.yLabel)
        ylabel(opt.yLabel, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    end
    if ~isempty(opt.title)
        title(opt.title, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
    end

    set(gca, 'FontSize', 18, 'FontName', 'Arial', 'TickDir', 'out', 'FontWeight', 'bold');
    ax = gca;
    ax.XAxis.FontSize = 18;
    ax.YAxis.FontSize = 18;
    ax.YAxis.FontWeight = 'bold';
    ax.LineWidth = 1.7;
    ax.Position = [0.18, 0.14, 0.78, 0.74];
    ax.XAxis.TickLength = [0.014 0.014];
    box off;

    % Legend for paired-line colour coding (shown if both classes present).
    legendHandles = gobjects(0);
    legendLabels  = {};
    if isgraphics(hIncrease)
        legendHandles(end+1) = hIncrease;
        legendLabels{end+1}  = 'Treatment > baseline';
    end
    if isgraphics(hDecrease)
        legendHandles(end+1) = hDecrease;
        legendLabels{end+1}  = 'Treatment \leq baseline';
    end
    if ~isempty(legendHandles)
        legend(legendHandles, legendLabels, 'Location', 'northwest', ...
            'Box', 'off', 'FontSize', 12);
    end

    if ~isempty(statsString)
        % Normalized axes coords so the stats box is independent of the
        % data range. Legend lives at 'northwest'; stats at top-right.
        text(0.985, 0.985, statsString, ...
            'Units',               'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment',   'top', ...
            'FontSize', 11, 'FontName', 'Arial', ...
            'BackgroundColor', [1 1 1], ...
            'EdgeColor', 'none', 'Margin', 4);
    end

    save_figure(fig, outFile);
    close(fig);
end

% -------------------------------------------------------------------------
function ymax = robust_ymax(y, mode)
% Pick a y-axis ceiling that is not dominated by a handful of outliers.
%   'auto'   : classic max-of-data
%   'robust' : max(97.5th percentile * 1.3, IQR upper fence, median * 4)
    switch lower(mode)
        case 'auto'
            ymax = max(y, [], 'omitnan');
        case 'robust'
            y = y(isfinite(y));
            if isempty(y)
                ymax = 1;
                return;
            end
            q1 = prctile(y, 25);
            q3 = prctile(y, 75);
            upperFence = q3 + 1.5 * (q3 - q1);
            p975 = prctile(y, 97.5) * 1.3;
            medFloor = max(0.001, median(y) * 4);
            ymax = max([upperFence, p975, medFloor]);
            dataMax = max(y);
            if ymax > dataMax
                ymax = dataMax;
            end
            if ymax <= 0
                ymax = max(y);
            end
        otherwise
            ymax = max(y, [], 'omitnan');
    end
    if ymax <= 0 || ~isfinite(ymax)
        ymax = 1;
    end
end
