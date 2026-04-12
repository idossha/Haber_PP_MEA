function plot_mean_iqr_marker(ax, x, yData, capHalfWidth, lineWidth, meanMarkerSize)
%PLOT_MEAN_IQR_MARKER Draw a vertical Q1-Q3 bar with caps and a mean dot.
%
%   PLOT_MEAN_IQR_MARKER(ax, x, yData, capHalfWidth, lineWidth, meanMarkerSize)
%   adds a black IQR marker to the supplied axes at horizontal position x.
%   The bar spans Q1-Q3 and the dot marks the mean. The same Tukey 1.5*IQR
%   inlier rule used by the boxplot is applied first so the bar matches the
%   boxplot whiskers visually.
%
% INPUTS:
%   ax              -  Target axes handle.
%   x               -  Horizontal position (axes data units).
%   yData           -  Vector of values for which to compute the marker.
%   capHalfWidth    -  Horizontal half-width of the Q1/Q3 caps (data units).
%   lineWidth       -  Line width for the IQR bar and caps.
%   meanMarkerSize  -  Marker size (points) for the mean dot.

    y = yData(:);
    q1f = prctile(y, 25);
    q3f = prctile(y, 75);
    iqr = q3f - q1f;
    if iqr > 0
        lo  = q1f - 1.5 * iqr;
        hi  = q3f + 1.5 * iqr;
        yIn = y(y >= lo & y <= hi);
    else
        yIn = y;
    end
    if isempty(yIn)
        yIn = y;
    end

    q1q3 = prctile(yIn, [25 75]);
    q1   = q1q3(1);
    q3   = q1q3(2);
    mu   = mean(yIn);

    hold(ax, 'on');
    h1 = plot(ax, [x x], [q1 q3], 'k-', 'LineWidth', lineWidth);
    h2 = plot(ax, [x - capHalfWidth, x + capHalfWidth], [q1 q1], 'k-', 'LineWidth', lineWidth);
    h3 = plot(ax, [x - capHalfWidth, x + capHalfWidth], [q3 q3], 'k-', 'LineWidth', lineWidth);
    h4 = plot(ax, x, mu, 'o', ...
        'Color', 'k', 'MarkerFaceColor', 'k', ...
        'MarkerSize', meanMarkerSize, 'LineWidth', 0.35);
    uistack([h1, h2, h3, h4], 'top');
end
