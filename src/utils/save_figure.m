function [pdfPath, pngPath] = save_figure(fig, outPath, varargin)
%SAVE_FIGURE Export a figure to a vector PDF (and a raster PNG shadow).
%
%   save_figure(fig, outPath) writes FIG as a **vector PDF** alongside a
%   raster PNG at the requested DPI (default 600). The PDF is the primary
%   artifact — it drops into Inkscape / Illustrator with every line,
%   patch, and text element still editable. The PNG is a raster shadow
%   for quick visual inspection and for tools that can't render PDFs
%   directly.
%
%   save_figure(fig, outPath, 'dpi', 600, 'format', {'pdf','png'})
%   overrides defaults. outPath may be passed with or without an
%   extension; both .pdf and .png companions are written next to it.
%
% INPUTS:
%   fig      -  Figure handle.
%   outPath  -  Output path. Extension is optional; the function strips
%               any extension and writes <base>.pdf and <base>.png.
%
% Name-value options:
%   'dpi'     -  Raster resolution for the PNG shadow and for any
%                rasterised elements inside the PDF (default 600).
%   'format'  -  Cell array of formats to write. Default {'pdf','png'}.
%                Pass {'pdf'} for PDF only, {'png'} for PNG only.
%
% OUTPUTS:
%   pdfPath  -  Absolute path to the written PDF (empty if not written).
%   pngPath  -  Absolute path to the written PNG (empty if not written).
%
% Notes:
%   - exportgraphics(..., 'ContentType','vector') keeps lines and patches
%     as vector geometry. Heatmaps / imagesc stay rasterised inside the
%     PDF at the requested DPI, which is correct behaviour.
%   - The 600-dpi PDF is the canonical deliverable; drop it into
%     Inkscape to recolour, relabel, or assemble into a composite figure.

    p = inputParser;
    addRequired(p, 'fig');
    addRequired(p, 'outPath');
    addParameter(p, 'dpi',    600, @(x) isscalar(x) && x >= 72);
    addParameter(p, 'format', {'pdf','png'}, @iscell);
    parse(p, fig, outPath, varargin{:});
    opt = p.Results;

    [parent, base, ~] = fileparts(outPath);
    if isempty(parent); parent = pwd; end
    if ~exist(parent, 'dir'); mkdir(parent); end

    pdfPath = '';
    pngPath = '';

    if any(strcmpi(opt.format, 'pdf'))
        pdfPath = fullfile(parent, [base '.pdf']);
        try
            exportgraphics(fig, pdfPath, ...
                'ContentType', 'vector', ...
                'Resolution',  opt.dpi, ...
                'BackgroundColor', 'white');
        catch ME
            warning('save_figure:PdfFallback', ...
                'Vector PDF failed (%s); falling back to print -dpdf.', ME.message);
            set(fig, 'PaperPositionMode', 'auto');
            print(fig, pdfPath, '-dpdf', sprintf('-r%d', opt.dpi));
        end
    end

    if any(strcmpi(opt.format, 'png'))
        pngPath = fullfile(parent, [base '.png']);
        try
            exportgraphics(fig, pngPath, 'Resolution', opt.dpi);
        catch ME
            warning('save_figure:PngFallback', ...
                'exportgraphics PNG failed (%s); falling back to print.', ME.message);
            set(fig, 'PaperPositionMode', 'auto');
            print(fig, pngPath, '-dpng', sprintf('-r%d', opt.dpi));
        end
    end
end
