function fig_filtered_traces(baselineName, treatmentName, channelIndices, varargin)
%FIG_FILTERED_TRACES Side-by-side filtered voltage traces (zoom + full).
%
%   FIG_FILTERED_TRACES(baselineName, treatmentName, channelIndices) loads
%   the requested baseline and treatment TDT blocks (paths resolved through
%   project_config()), filters every requested channel using the canonical
%   notch + bandpass stack, and renders a four-panel publication figure:
%
%       [ Baseline zoom | Baseline full | Treatment zoom | Treatment full ]
%
%   Each row is a channel, scaled symmetrically around its own per-channel
%   peak amplitude. Per-channel voltage scale bars sit in the right margin
%   and time scale bars (10 s zoom, 60 s full) sit below the panels. The
%   figure is written to output/shared/.
%
%   FIG_FILTERED_TRACES(..., 'labels', {'Baseline','DOI'}, ...
%                            'zoomDuration', 30)
%
% INPUTS:
%   baselineName    -  Dataset folder name for the baseline block.
%   treatmentName   -  Dataset folder name for the treatment block.
%   channelIndices  -  Vector of channel numbers to plot (top to bottom).
%
% Name-value options:
%   'labels'        -  {leftLabel, rightLabel}, default {'Baseline','DOI'}
%   'zoomDuration'  -  seconds shown in the zoom panels, default 30
%   'outFile'       -  override the default output path
%
% OUTPUTS:
%   PNG file written to output/shared/ (or 'outFile' if supplied).

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'baselineName',   @(s) ischar(s) || isstring(s));
    addRequired(p,  'treatmentName',  @(s) ischar(s) || isstring(s));
    addRequired(p,  'channelIndices', @(v) isnumeric(v) && ~isempty(v));
    addParameter(p, 'labels',         {'Baseline','DOI'});
    addParameter(p, 'zoomDuration',   30);
    addParameter(p, 'outFile',        '');
    parse(p, baselineName, treatmentName, channelIndices, varargin{:});
    opt = p.Results;

    addpath(genpath(cfg.paths.tdt_sdk));

    controlPath   = dataset_path(char(opt.baselineName),  cfg);
    treatmentPath = dataset_path(char(opt.treatmentName), cfg);

    nCh = numel(opt.channelIndices);
    fprintf('Loading %d channels for filtered-trace figure...\n', nCh);

    controlSignals   = cell(nCh, 1);
    treatmentSignals = cell(nCh, 1);
    fsControl   = NaN; fsTreatment = NaN;
    timeControl = []; timeTreatment = [];

    for i = 1:nCh
        ch = opt.channelIndices(i);
        [xC, fsC] = load_and_filter(controlPath,   ch, cfg);
        [xT, fsT] = load_and_filter(treatmentPath, ch, cfg);
        controlSignals{i}   = xC * 1e6;       % volts -> microvolts
        treatmentSignals{i} = xT * 1e6;
        if i == 1
            fsControl   = fsC;
            fsTreatment = fsT;
            timeControl   = (0:numel(xC)-1) / fsC;
            timeTreatment = (0:numel(xT)-1) / fsT;
        end
    end

    if isempty(opt.outFile)
        sharedDir = output_path(cfg, '', 'shared', '');
        if ~exist(sharedDir, 'dir')
            mkdir(sharedDir);
        end
        outFile = fullfile(sharedDir, ...
            sprintf('filtered_traces_%s_vs_%s.png', ...
                strrep(char(opt.baselineName),  '-', '_'), ...
                strrep(char(opt.treatmentName), '-', '_')));
    else
        outFile = opt.outFile;
    end

    render_filtered_traces(controlSignals, treatmentSignals, ...
        timeControl, timeTreatment, opt.zoomDuration, opt.labels, outFile);

    fprintf('fig_filtered_traces: saved %s\n', outFile);
end

% =========================================================================
function render_filtered_traces(cSig, tSig, tC, tT, zoomDuration, ~, outFile)
% Layout, scale bars, and styling are kept identical to the
% proof-of-concept figures_compare_filtered_signals.m output.

    nChannels  = numel(cSig);
    scaleBarZoom    = 10;
    scaleBarFull    = 60;
    xScaleBarLineWidth = 4.5;
    xScaleBarCapWidth  = 3.5;
    scaleBarFontSize   = 21;
    figPosWidth        = 2400;
    figPaperWidth      = 22;
    midGapBaselineDOI  = 0.035;

    fig = figure('Visible', 'off', 'Position', [100 100 figPosWidth 1000], 'Color', 'w');
    set(fig, 'Units', 'inches', 'PaperUnits', 'inches', ...
        'PaperSize', [figPaperWidth 10], 'PaperPosition', [0 0 figPaperWidth 10]);

    boxHalfHeight = 250;
    channelGap    = -95;
    channelOffset = 2 * boxHalfHeight + channelGap;

    yPositions = zeros(nChannels, 1);
    for i = 1:nChannels
        yPositions(i) = (nChannels - i) * channelOffset;
    end
    bottomBoxBottom = min(yPositions) - boxHalfHeight;
    scaleBarY      = bottomBoxBottom - 15;
    xScaleLabelOff = 12;
    yLimMin = scaleBarY - xScaleLabelOff - 22;
    yLimMax = max(yPositions) + boxHalfHeight + 25;

    maxTC = max(tC); maxTT = max(tT);
    idxCZ = tC <= min(zoomDuration, maxTC);
    idxTZ = tT <= min(zoomDuration, maxTT);

    maxAbs = zeros(nChannels, 1);
    for i = 1:nChannels
        s1 = cSig{i}(idxCZ);
        s2 = cSig{i};
        s3 = tSig{i}(idxTZ);
        s4 = tSig{i};
        A_i = max([max(abs(s1)), max(abs(s2)), max(abs(s3)), max(abs(s4))]);
        if A_i < 1, A_i = 20; end
        maxAbs(i) = A_i;
    end

    ax_bz = axes('Position', [0.12, 0.06, 0.095, 0.88]);
    draw_panel(cSig, tC, idxCZ, yPositions, maxAbs, boxHalfHeight, true);
    finalize_panel(ax_bz, [0, min(zoomDuration, maxTC)], [yLimMin yLimMax]);
    draw_time_scale_bar(scaleBarZoom, scaleBarY, xScaleLabelOff, '10 s', ...
        scaleBarFontSize, xScaleBarLineWidth, xScaleBarCapWidth);

    ax_bf = axes('Position', [0.225, 0.06, 0.29, 0.88]);
    draw_panel(cSig, tC, true(size(tC)), yPositions, maxAbs, boxHalfHeight, false);
    finalize_panel(ax_bf, [0, maxTC], [yLimMin yLimMax]);
    draw_time_scale_bar(scaleBarFull, scaleBarY, xScaleLabelOff, '60 s', ...
        scaleBarFontSize, xScaleBarLineWidth, xScaleBarCapWidth);

    ax_tz = axes('Position', [0.225 + 0.29 + midGapBaselineDOI, 0.06, 0.095, 0.88]);
    draw_panel(tSig, tT, idxTZ, yPositions, maxAbs, boxHalfHeight, true);
    finalize_panel(ax_tz, [0, min(zoomDuration, maxTT)], [yLimMin yLimMax]);
    draw_time_scale_bar(scaleBarZoom, scaleBarY, xScaleLabelOff, '10 s', ...
        scaleBarFontSize, xScaleBarLineWidth, xScaleBarCapWidth);

    ax_tf = axes('Position', [0.225 + 0.29 + midGapBaselineDOI + 0.095 + 0.015, 0.06, 0.27, 0.88]);
    draw_panel(tSig, tT, true(size(tT)), yPositions, maxAbs, boxHalfHeight, false);
    finalize_panel(ax_tf, [0, maxTT], [yLimMin yLimMax]);
    draw_time_scale_bar(scaleBarFull, scaleBarY, xScaleLabelOff, '60 s', ...
        scaleBarFontSize, xScaleBarLineWidth, xScaleBarCapWidth);

    % Voltage scale bars in the right margin.
    axPos    = get(ax_bz, 'Position');
    axRight  = get(ax_tf, 'Position');
    rightEdgeTraces = axRight(1) + axRight(3);
    marginBar       = 0.012;
    barCenterX      = rightEdgeTraces + marginBar;
    barHalfLen      = 0.006;
    voltLabelNormHeight = 0.032;
    textGapFromBar      = 0.004;
    voltLabelNormWidth  = min(0.055, 0.998 - (barCenterX + textGapFromBar));
    for i = 1:nChannels
        ampRound = round(2 * maxAbs(i) / 20) * 20;
        if ampRound < 20, ampRound = 20; end
        yn = axPos(2) + axPos(4) * (yPositions(i) - yLimMin) / (yLimMax - yLimMin);
        annotation(fig, 'line', [barCenterX, barCenterX], ...
            [yn - barHalfLen, yn + barHalfLen], 'Color', 'k', 'LineWidth', 2);
        annotation(fig, 'textbox', ...
            [barCenterX + textGapFromBar, yn - voltLabelNormHeight/2, ...
             voltLabelNormWidth, voltLabelNormHeight], ...
            'String', sprintf('%d µV', ampRound), ...
            'FontSize', scaleBarFontSize, 'FontName', 'Arial', 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'VerticalAlignment', 'middle', ...
            'Interpreter', 'none', 'FitBoxToText', 'on');
    end

    save_figure(fig, outFile);
    close(fig);
end

% -------------------------------------------------------------------------
function draw_panel(signals, timeVec, idxKeep, yPositions, maxAbs, boxHalfHeight, isZoom)
    hold on;
    nChannels = numel(signals);
    targetSamples = 15000;
    if ~isZoom
        targetSamples = 50000;
    end
    for i = 1:nChannels
        sig = signals{i};
        tk  = timeVec(idxKeep);
        sk  = sig(idxKeep);
        ds  = max(1, floor(numel(sk) / targetSamples));
        tp  = tk(1:ds:end);
        sp  = sk(1:ds:end);
        A_i = maxAbs(i);
        yp  = yPositions(i) + (sp / A_i) * boxHalfHeight;
        plot(tp, yp, 'Color', [0 0 0], 'LineWidth', 0.65);
    end
end

% -------------------------------------------------------------------------
function finalize_panel(ax, xLimVec, yLimVec)
    xlim(ax, xLimVec); ylim(ax, yLimVec);
    set(ax, 'XTick', [], 'XTickLabel', [], 'YTick', [], 'YTickLabel', []);
    set(ax, 'FontSize', 10, 'FontName', 'Arial', 'TickDir', 'out', ...
        'XColor', 'none', 'YColor', 'none');
    box(ax, 'off'); grid(ax, 'off');
end

% -------------------------------------------------------------------------
function draw_time_scale_bar(barLen, yBar, labelOff, labelText, fontSize, lwBar, lwCap)
    plot([0, barLen], [yBar, yBar], 'k-', 'LineWidth', lwBar);
    plot([0, 0], [yBar - 6, yBar + 6], 'k-', 'LineWidth', lwCap);
    plot([barLen, barLen], [yBar - 6, yBar + 6], 'k-', 'LineWidth', lwCap);
    text(barLen/2, yBar - labelOff, labelText, 'FontSize', fontSize, ...
        'FontName', 'Arial', 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
