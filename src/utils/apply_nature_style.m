function apply_nature_style(fig)
% APPLY_NATURE_STYLE  Enforce Nature/NPP figure specifications.
%
%   apply_nature_style(fig) sets font, line, tick, and annotation properties
%   on figure FIG and ALL child objects to match Nature portfolio artwork
%   requirements:
%
%     Tick labels     6 pt Arial
%     Axis labels     7 pt Arial  (via LabelFontSizeMultiplier)
%     Titles          7 pt Arial bold
%     All other text  capped at 7 pt, floored at 5 pt
%     Axes frame      0.6 pt
%     Data lines      floored at 0.25 pt, capped at 1.5 pt
%     Annotations     line weight floored at 0.25 pt, capped at 1.5 pt
%     Ticks           outward, no box, no grid
%
%   This is the SINGLE source of truth for styling.  Figure scripts should
%   not set font sizes or axis line widths manually — call this function
%   once, immediately before save_figure.
%
%   Nature specs (https://research-figure-guide.nature.com/figures/):
%     - Font: Arial or Helvetica, 5-7 pt body, 8 pt bold lowercase panel labels
%     - Minimum line weight: 0.25 pt
%     - Axes: 0.6 pt
%     - Color: RGB
%     - Panel labels: 8 pt bold lowercase a, b, c
%
% See also: CREATE_PANEL_FIGURE, SAVE_FIGURE.

    if nargin < 1; fig = gcf; end

    % --- Constants (Nature portfolio artwork guide) ---
    FONT_NAME       = 'Arial';
    TICK_SIZE       = 6;            % pt — tick labels
    LABEL_SIZE      = 7;            % pt — axis labels, titles
    MIN_TEXT_SIZE   = 5;            % pt — Nature minimum readable
    MAX_TEXT_SIZE   = 7;            % pt — ceiling for body text (excl. 8pt panel labels)
    AXES_LINE_WIDTH = 0.6;          % pt — axes frame and ticks
    MIN_LINE_WIDTH  = 0.25;         % pt — Nature minimum line weight
    MAX_DATA_LINE   = 1.5;          % pt — clamp for oversized data lines
    MAX_ANNO_LINE   = 1.5;          % pt — clamp for annotation lines (scale bars etc.)
    TICK_LENGTH     = [0.015 0.025];

    % --- Figure-level defaults ---
    set(fig, 'Color', 'w', 'Renderer', 'painters');

    % --- Axes ---
    allAxes = findall(fig, 'Type', 'axes');
    for i = 1:numel(allAxes)
        ax = allAxes(i);

        % Font: tick labels at TICK_SIZE; labels and titles at LABEL_SIZE
        set(ax, 'FontName', FONT_NAME, 'FontSize', TICK_SIZE);
        set(ax, 'LabelFontSizeMultiplier', LABEL_SIZE / TICK_SIZE);
        set(ax, 'TitleFontSizeMultiplier', LABEL_SIZE / TICK_SIZE);
        set(ax, 'TitleFontWeight', 'bold');

        % Frame and ticks
        set(ax, 'LineWidth', AXES_LINE_WIDTH);
        set(ax, 'TickLength', TICK_LENGTH);
        set(ax, 'TickDir', 'out');
        set(ax, 'Box', 'off');
        set(ax, 'XGrid', 'off', 'YGrid', 'off');

        % Enforce line weight bounds on all lines within this axes
        lines = findall(ax, 'Type', 'line');
        for j = 1:numel(lines)
            lw = get(lines(j), 'LineWidth');
            if lw < MIN_LINE_WIDTH
                set(lines(j), 'LineWidth', MIN_LINE_WIDTH);
            elseif lw > MAX_DATA_LINE
                set(lines(j), 'LineWidth', MAX_DATA_LINE);
            end
        end
    end

    % --- All line-bearing objects (including annotations) ---
    allLines = findall(fig, 'Type', 'line');
    for i = 1:numel(allLines)
        lw = get(allLines(i), 'LineWidth');
        if lw < MIN_LINE_WIDTH
            set(allLines(i), 'LineWidth', MIN_LINE_WIDTH);
        elseif lw > MAX_ANNO_LINE
            set(allLines(i), 'LineWidth', MAX_ANNO_LINE);
        end
    end

    % Also enforce on annotation lines (they use a different object type)
    allAnno = findall(fig, 'Type', 'annotationline', '-or', 'Type', 'line');
    for i = 1:numel(allAnno)
        try
            lw = get(allAnno(i), 'LineWidth');
            if lw < MIN_LINE_WIDTH
                set(allAnno(i), 'LineWidth', MIN_LINE_WIDTH);
            elseif lw > MAX_ANNO_LINE
                set(allAnno(i), 'LineWidth', MAX_ANNO_LINE);
            end
        catch
        end
    end

    % --- Patch objects (scatter fills, rectangles, etc.) ---
    allPatches = findall(fig, 'Type', 'patch');
    for i = 1:numel(allPatches)
        try
            lw = get(allPatches(i), 'LineWidth');
            if lw > 0 && lw < MIN_LINE_WIDTH
                set(allPatches(i), 'LineWidth', MIN_LINE_WIDTH);
            end
        catch
        end
    end

    % --- All text-bearing objects ---
    allObj = findall(fig, '-property', 'FontSize');
    for i = 1:numel(allObj)
        obj = allObj(i);
        if isa(obj, 'matlab.graphics.axis.Axes'); continue; end
        if isa(obj, 'matlab.ui.Figure');          continue; end
        try
            set(obj, 'FontName', FONT_NAME);
            fs = get(obj, 'FontSize');
            if fs > MAX_TEXT_SIZE
                set(obj, 'FontSize', MAX_TEXT_SIZE);
            elseif fs < MIN_TEXT_SIZE && fs > 0
                set(obj, 'FontSize', MIN_TEXT_SIZE);
            end
        catch
        end
    end
end
