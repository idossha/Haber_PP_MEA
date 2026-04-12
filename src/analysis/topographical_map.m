function fig = topographical_map(channelMetric, varargin)
%TOPOGRAPHICAL_MAP Render a per-electrode metric over the 60-channel MEA grid.
%
%   fig = TOPOGRAPHICAL_MAP(channelMetric) draws an 8x8 heatmap of the
%   supplied per-channel metric using the standard 60-electrode layout
%   (corners blanked). The order of channelMetric is the recording channel
%   index 1..nCh; the channel-to-grid mapping is read from the shared
%   mea60_layout() helper (src/utils/mea60_layout.m), which is also used
%   by plot_network_on_mea.m so every spatial plot in the codebase shares
%   the same geometry.
%
% INPUTS:
%   channelMetric  -  Vector of length nCh (default expectation: 60). Use
%                     NaN for channels you want masked white.
%
% Name-value options:
%   'title'      -  Figure title (default '').
%   'colormap'   -  Colormap name or matrix (default 'parula').
%   'climits'    -  [lo hi] color limits (default = data range, NaNs ignored).
%   'cbarLabel'  -  Colorbar label (default '').
%   'outFile'    -  If non-empty, save the figure as PNG.
%
% OUTPUTS:
%   fig  -  Figure handle.
%
% Notes:
%   - The 60MEA100/10iR layout is an 8x8 grid with the four corners absent
%     and one reference electrode replaced. We blank cells corresponding to
%     missing positions with NaN before plotting.
%   - For studies that use the full 64-channel TDT acquisition you can pass
%     a length-64 vector; only positions present in the layout map are
%     drawn, the remainder are masked.

    p = inputParser;
    addRequired(p,  'channelMetric', @(v) isnumeric(v) && isvector(v));
    addParameter(p, 'title',     '');
    addParameter(p, 'colormap',  'parula');
    addParameter(p, 'climits',   []);
    addParameter(p, 'cbarLabel', '');
    addParameter(p, 'outFile',   '');
    parse(p, channelMetric, varargin{:});
    opt = p.Results;

    layout = mea60_layout();   % shared geometry (src/utils/mea60_layout.m)

    grid = nan(layout.gridSize);
    for k = 1:layout.nCh
        ch = layout.channelList(k);
        if ch > numel(channelMetric); continue; end
        grid(layout.gridRow(k), layout.gridCol(k)) = channelMetric(ch);
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 600 600]);
    imAlpha = ~isnan(grid);
    imagesc(grid, 'AlphaData', imAlpha);
    axis image off;
    colormap(gca, opt.colormap);
    if isempty(opt.climits)
        finiteVals = grid(~isnan(grid));
        if ~isempty(finiteVals)
            clim([min(finiteVals), max(finiteVals)]);
        end
    else
        clim(opt.climits);
    end
    drawnow limitrate nocallbacks;
    cb = colorbar;
    if ~isempty(opt.cbarLabel)
        cb.Label.String     = opt.cbarLabel;
        cb.Label.FontSize   = 14;
        cb.Label.FontName   = 'Arial';
        cb.Label.FontWeight = 'bold';
    end
    if ~isempty(opt.title)
        title(opt.title, 'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');
    end

    if ~isempty(opt.outFile)
        save_figure(fig, opt.outFile);
    end
end
