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
%     All other text  capped at 7 pt Arial (legends, annotations, etc.)
%     Axes frame      0.6 pt
%     Data lines      capped at 1.0 pt (lines thicker than 1.5 pt are reduced)
%     Ticks           outward, no box, no grid
%
%   This is the SINGLE source of truth for styling.  Figure scripts should
%   not set font sizes or axis line widths manually — call this function
%   once, immediately before save_figure.
%
% See also: CREATE_PANEL_FIGURE, SAVE_FIGURE.

    if nargin < 1; fig = gcf; end

    % --- Constants (Nature/NPP artwork guide) ---
    FONT_NAME       = 'Arial';
    TICK_SIZE       = 6;            % pt — tick labels
    LABEL_SIZE      = 7;            % pt — axis labels, titles
    MAX_TEXT_SIZE   = 7;            % pt — ceiling for all non-tick text
    AXES_LINE_WIDTH = 0.6;          % pt — axes frame
    MAX_DATA_LINE   = 1.0;          % pt — clamp for oversized data lines
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

        % Clamp oversized data lines within this axes
        lines = findall(ax, 'Type', 'line');
        for j = 1:numel(lines)
            if get(lines(j), 'LineWidth') > 1.5
                set(lines(j), 'LineWidth', MAX_DATA_LINE);
            end
        end
    end

    % --- All text-bearing objects ---
    % Using '-property' finds EVERY object that has a FontSize property,
    % including annotation textboxes, colorbar labels, and legend entries
    % that findall(..., 'Type', 'text') misses.
    allObj = findall(fig, '-property', 'FontSize');
    for i = 1:numel(allObj)
        obj = allObj(i);
        % Skip axes and figure (already styled above)
        if isa(obj, 'matlab.graphics.axis.Axes'); continue; end
        if isa(obj, 'matlab.ui.Figure');          continue; end
        try
            set(obj, 'FontName', FONT_NAME);
            if get(obj, 'FontSize') > MAX_TEXT_SIZE
                set(obj, 'FontSize', MAX_TEXT_SIZE);
            end
        catch
            % Some internal graphics objects do not support direct setting
        end
    end
end
