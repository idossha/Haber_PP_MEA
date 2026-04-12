function plot_paired_boxplot(yBaseline, yTreatment, colors, labels, outFile, varargin)
%PLOT_PAIRED_BOXPLOT Two-group boxplot with mean overlay (Tukey inlier rule).
%
%   PLOT_PAIRED_BOXPLOT(yBaseline, yTreatment, colors, labels, outFile)
%   renders the canonical paired boxplot used by fig_spike_rate /
%   fig_burst_rate. Boxes use Tukey 1.5*IQR whiskers; the overlaid mean is
%   computed on inliers only.
%
%   PLOT_PAIRED_BOXPLOT(..., 'Name', value, ...) with any of:
%     'title'       figure title
%     'yLabel'      y-axis label (e.g. 'Spike rate (spikes/min)')
%     'yLimitMode'  'auto' | 'robust' (default 'robust')
%
% INPUTS:
%   yBaseline   -  Column vector of baseline values.
%   yTreatment  -  Matching column vector of treatment values.
%   colors      -  Struct from paired_plot_colors().
%   labels      -  Struct with .baseline / .treatment label strings.
%   outFile     -  Absolute path for the saved PNG.

    p = inputParser;
    addRequired(p,  'yBaseline');
    addRequired(p,  'yTreatment');
    addRequired(p,  'colors');
    addRequired(p,  'labels');
    addRequired(p,  'outFile');
    addParameter(p, 'title',      '');
    addParameter(p, 'yLabel',     '');
    addParameter(p, 'yLimitMode', 'robust', @(s) any(strcmpi(s, {'auto','robust'})));
    addParameter(p, 'yScale',     'linear', @(s) any(strcmpi(s, {'linear','log'})));
    addParameter(p, 'yFloor',     0.3,      @(x) isscalar(x) && x > 0);
    parse(p, yBaseline, yTreatment, colors, labels, outFile, varargin{:});
    opt = p.Results;

    nCh = numel(yBaseline);
    if strcmpi(opt.yScale, 'log')
        yBase  = max(yBaseline,  opt.yFloor);
        yTreat = max(yTreatment, opt.yFloor);
    else
        yBase  = yBaseline;
        yTreat = yTreatment;
    end
    y = [yBase; yTreat];
    g = [ones(nCh, 1); 2 * ones(nCh, 1)];

    fig = figure('Visible', 'off', 'Position', [100, 100, 720, 720], 'Color', 'w');
    set(fig, 'Units', 'inches', 'PaperUnits', 'inches', ...
        'PaperSize', [7.2 7.2], 'PaperPosition', [0 0 7.2 7.2]);
    hold on;

    boxplot(y, g, ...
        'Positions', [1, 1.26], ...
        'Labels',    {labels.baseline, labels.treatment}, ...
        'Widths',    0.14, ...
        'Symbol',    '');

    hBox = findobj(gca, 'Tag', 'Box');
    for j = 1:length(hBox)
        xd = get(hBox(j), 'XData');
        yd = get(hBox(j), 'YData');
        if j == 1
            hp = patch(xd, yd, colors.treatmentBox, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
            uistack(hp, 'bottom');
            set(hBox(j), 'Color', colors.treatmentOutline, 'LineWidth', 1.2);
        else
            hp = patch(xd, yd, colors.baselineBox, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
            uistack(hp, 'bottom');
            set(hBox(j), 'Color', colors.baselineOutline, 'LineWidth', 1.2);
        end
    end
    hMed = findobj(gca, 'Tag', 'Median');
    set(hMed, 'Color', [0 0 0], 'LineWidth', 1.5);
    hWhisker = [findobj(gca, 'Tag', 'Upper Whisker'); findobj(gca, 'Tag', 'Lower Whisker')];
    if ~isempty(hWhisker)
        set(hWhisker, 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 1.2);
    end

    hCap = [findobj(gca, 'Tag', 'Upper Adjacent Value'); findobj(gca, 'Tag', 'Lower Adjacent Value')];
    for k = 1:length(hCap)
        xd = get(hCap(k), 'XData');
        xm = (xd(1) + xd(2)) / 2;
        d  = (xd(2) - xd(1)) * 0.5;
        set(hCap(k), 'XData', [xm - d/2, xm + d/2], 'LineWidth', 1.2);
    end

    % Inlier mean overlay (Tukey 1.5*IQR rule).
    meanBaseline  = inlier_mean(yBase);
    meanTreatment = inlier_mean(yTreat);
    xd1 = get(hBox(1), 'XData'); x1L = min(xd1); x1R = max(xd1);
    xd2 = get(hBox(2), 'XData'); x2L = min(xd2); x2R = max(xd2);
    hMean1 = plot([x2L, x2R], [meanBaseline,  meanBaseline ], '--', ...
        'Color', colors.baselineOutline,  'LineWidth', 1.2);
    hMean2 = plot([x1L, x1R], [meanTreatment, meanTreatment], '--', ...
        'Color', colors.treatmentOutline, 'LineWidth', 1.2);
    uistack([hMean1, hMean2], 'top');

    xlim([0.82, 1.38]);
    xticks([1, 1.26]);
    xticklabels({labels.baseline, labels.treatment});

    % --- Robust y-axis --------------------------------------------------
    ymax = robust_ymax([yBase; yTreat], opt.yLimitMode);
    if strcmpi(opt.yScale, 'log')
        set(gca, 'YScale', 'log');
        yTop = ymax * 1.20;
        if yTop <= opt.yFloor
            yTop = opt.yFloor * 10;
        end
        ylim([opt.yFloor, yTop]);
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
    ax.Position = [0.18, 0.14, 0.78, 0.76];
    ax.XAxis.TickLength = [0.014 0.014];
    box off;

    save_figure(fig, outFile);
    close(fig);
end

% -------------------------------------------------------------------------
function mu = inlier_mean(y)
    q1 = prctile(y, 25);
    q3 = prctile(y, 75);
    iqr = q3 - q1;
    inRange = y >= (q1 - 1.5 * iqr) & y <= (q3 + 1.5 * iqr);
    if any(inRange)
        mu = mean(y(inRange));
    else
        mu = mean(y);
    end
end

% -------------------------------------------------------------------------
function ymax = robust_ymax(y, mode)
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
