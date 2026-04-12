function plot_pct_change_violin(pctData, colors, xLabelText, yAxisLabel, pctYlim, outFile)
%PLOT_PCT_CHANGE_VIOLIN Vertical violin + boxplot of percent-change values.
%
%   PLOT_PCT_CHANGE_VIOLIN(pctData, colors, xLabelText, yAxisLabel, pctYlim, outFile)
%   renders the canonical pct-change violin used by
%   fig_spike_pct_change_violin and fig_burst_pct_change_violin. The y-axis
%   floor is fixed at -100%; the upper limit is rounded up to the nearest
%   100% from the Tukey 1.5*IQR upper fence (with data safeguard) when
%   pctYlim is empty, otherwise pctYlim(2) is rounded up to the next 100%.
%   Box uses 1.5*IQR whiskers (no symbols), median is dashed, and the mean
%   line uses the Tukey-inlier rule.
%
% INPUTS:
%   pctData     -  Column vector of percent-change observations.
%   colors      -  Struct from paired_plot_colors() (treatment + outline).
%   xLabelText  -  Single x-axis label (e.g. 'DOI', 'Ketanserin').
%   yAxisLabel  -  y-axis label (e.g. 'Percent change in spike rate (%)').
%   pctYlim     -  Optional [yLo yHi]; pass [] for auto top.
%   outFile     -  Absolute path for the saved PNG.

    if isempty(pctData)
        error('plot_pct_change_violin:NoData', 'pctData is empty.');
    end

    histFaceAlpha     = 0.75;
    violinFaceColor   = 0.35 * colors.treatment       + 0.65 * [1 1 1];
    violinEdgeColor   = 0.55 * colors.treatmentOutline + 0.45 * [1 1 1];
    violinHalfWidth   = 0.32;
    violinYGridN      = 180;
    violinFaceAlpha   = 0.45;
    violinBottomGapFr = 0.001;
    xCenter           = 1;
    boxWidth          = 0.16;

    q1       = prctile(pctData, 25);
    q3       = prctile(pctData, 75);
    iqrPct   = q3 - q1;
    upper    = q3 + 1.5 * iqrPct;

    if isempty(pctYlim)
        yRound = ceil(upper / 100) * 100;
        maxPct = max(pctData);
        if maxPct > yRound
            yLimTop = ceil(maxPct / 100) * 100;
        else
            yLimTop = yRound;
        end
        labelHeadroom = max(15, 0.03 * (yLimTop + 100));
        yLimHi = ceil((yLimTop + labelHeadroom) / 100) * 100;
        yLimUse = [-100, yLimHi];
    else
        yLimUse    = pctYlim(:)';
        yLimUse(1) = -100;
        yLimUse(2) = ceil(yLimUse(2) / 100) * 100;
    end

    medPct  = median(pctData, 'omitnan');
    inRange = pctData >= (q1 - 1.5 * iqrPct) & pctData <= (q3 + 1.5 * iqrPct);
    if any(inRange)
        meanPctInlier = mean(pctData(inRange));
    else
        meanPctInlier = mean(pctData);
    end

    fprintf('  median %.2f%% | IQR [%.2f, %.2f]%% | mean(inlier) %.2f%% | n=%d\n', ...
        medPct, q1, q3, meanPctInlier, numel(pctData));

    fig = figure('Visible', 'off', 'Position', [100 100 380 620], 'Color', 'w');
    set(fig, 'Units', 'inches', 'PaperUnits', 'inches', ...
        'PaperSize', [3.8 6.2], 'PaperPosition', [0 0 3.8 6.2]);
    hold on;
    ylim(yLimUse);

    spanY    = yLimUse(2) - yLimUse(1);
    violinYLo = yLimUse(1) + max(0, violinBottomGapFr * spanY);
    yGrid    = linspace(violinYLo, yLimUse(2), violinYGridN)';
    try
        [f, ~] = ksdensity(pctData, yGrid, 'BoundaryCorrection', 'reflection');
    catch
        [f, ~] = ksdensity(pctData, yGrid);
    end

    hViolin = [];
    fMax = max(f);
    if fMax > 0
        fNorm = (f / fMax) * violinHalfWidth;
        xL = xCenter - fNorm(:);
        xR = xCenter + fNorm(:);
        hViolin = patch([xL; flipud(xR)], [yGrid; flipud(yGrid)], ...
            violinFaceColor, 'FaceAlpha', violinFaceAlpha, 'EdgeColor', 'none');
        uistack(hViolin, 'bottom');
        plot(xL, yGrid, 'Color', violinEdgeColor, 'LineWidth', 1);
        plot(xR, yGrid, 'Color', violinEdgeColor, 'LineWidth', 1);
        uistack(hViolin, 'bottom');
    end

    yline(0, '-', 'Color', [0 0 0], 'LineWidth', 1.0);

    boxplot(pctData, 'Positions', xCenter, 'Widths', boxWidth, ...
        'Symbol', '', 'Colors', [0 0 0]);
    ylim(yLimUse);

    hBox = findobj(gca, 'Tag', 'Box');
    for j = 1:numel(hBox)
        xd = get(hBox(j), 'XData');
        yd = get(hBox(j), 'YData');
        hpBox = patch(xd, yd, colors.treatment, 'FaceAlpha', histFaceAlpha, 'EdgeColor', 'none');
        uistack(hpBox, 'bottom');
        if ~isempty(hViolin)
            uistack(hViolin, 'bottom');
        end
        set(hBox(j), 'Color', colors.treatmentOutline, 'LineWidth', 1.2);
    end
    hMed = findobj(gca, 'Tag', 'Median');
    if ~isempty(hMed)
        set(hMed, 'Color', violinEdgeColor, 'LineWidth', 2, 'LineStyle', '--');
    end
    hWh = [findobj(gca, 'Tag', 'Upper Whisker'); findobj(gca, 'Tag', 'Lower Whisker')];
    if ~isempty(hWh)
        set(hWh, 'LineWidth', 1.1, 'Color', [0 0 0], 'LineStyle', '-');
    end
    hAdj = [findobj(gca, 'Tag', 'Upper Adjacent Value'); findobj(gca, 'Tag', 'Lower Adjacent Value')];
    for k = 1:numel(hAdj)
        xd = get(hAdj(k), 'XData');
        xm = (xd(1) + xd(2)) / 2;
        d  = (xd(2) - xd(1)) * 0.5;
        set(hAdj(k), 'XData', [xm - d/2, xm + d/2], 'LineWidth', 1.2, ...
            'Color', [0 0 0], 'LineStyle', '-');
    end

    if ~isempty(hBox)
        xdBox = get(hBox(1), 'XData');
        xMeanLeft  = min(xdBox);
        xMeanRight = max(xdBox);
        hMeanLn = plot([xMeanLeft, xMeanRight], [meanPctInlier, meanPctInlier], '-', ...
            'Color', [0 0 0], 'LineWidth', 1.2);
        uistack(hMeanLn, 'top');
    end

    xlim([0.62, 1.38]);
    xticks(xCenter); xticklabels({xLabelText});
    ylabel(yAxisLabel, 'FontSize', 18, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 18, 'FontName', 'Arial', 'TickDir', 'out', ...
        'FontWeight', 'bold', 'YGrid', 'off', 'XGrid', 'off');
    ax = gca;
    ax.XAxis.FontSize = 18;
    ax.YAxis.FontSize = 18;
    ax.YAxis.FontWeight = 'bold';
    ax.LineWidth = 1.7;
    ax.Position = [0.20, 0.10, 0.72, 0.74];
    ax.XAxis.TickLength = [0.014 0.014];
    ax.YAxis.TickLength = [0.014 0.014];
    box off;

    % Re-draw bottom axis line on top.
    hXAxis = plot([0.62, 1.38], [yLimUse(1), yLimUse(1)], '-', ...
        'Color', [0 0 0], 'LineWidth', 1.7);
    uistack(hXAxis, 'top');

    yTickStep = nice_tick_step(yLimUse(2) - yLimUse(1), 9);
    ytVals = yLimUse(1):yTickStep:yLimUse(2);
    if isempty(ytVals) || ytVals(end) < yLimUse(2)
        ytVals = [ytVals, yLimUse(2)];
    end
    ytVals = unique([ytVals(:); 0]);
    ytVals = sort(ytVals);
    yticks(ytVals);

    save_figure(fig, outFile);
    close(fig);
end
