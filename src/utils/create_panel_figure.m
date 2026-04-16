function fig = create_panel_figure(widthCm, heightCm)
%CREATE_PANEL_FIGURE Create a publication figure at exact physical dimensions.
%
%   fig = CREATE_PANEL_FIGURE(widthCm, heightCm) returns an invisible figure
%   whose screen aspect ratio matches the paper output size.  This prevents
%   the normalised-coordinate distortion that occurs when the screen figure
%   has a different aspect ratio to the PaperSize.
%
%   Usage pattern:
%       fig = create_panel_figure(5.8, 4.4);
%       ax  = axes(fig);
%       ... plotting code ...
%       apply_nature_style(fig);
%       save_figure(fig, outPath);
%
% See also: APPLY_NATURE_STYLE, SAVE_FIGURE.

    PX_PER_CM = 120;   % screen magnification (ensures legible layout)

    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100 100 round(widthCm * PX_PER_CM) ...
                             round(heightCm * PX_PER_CM)]);
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperSize', [widthCm heightCm]);
    set(fig, 'PaperPosition', [0 0 widthCm heightCm]);
    set(fig, 'Renderer', 'painters');
end
