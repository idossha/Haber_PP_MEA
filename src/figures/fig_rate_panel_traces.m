function fig_rate_panel_traces(study, varargin)
%FIG_RATE_PANEL_TRACES Filtered voltage traces for baseline vs treatment.
%
%   FIG_RATE_PANEL_TRACES(study) picks one representative pair from
%   STUDY ('doi' or 'ket'), loads and filters 10 SNR+story-ranked
%   channels from the baseline and treatment datasets, and draws a
%   panel A figure with four sub-columns per row:
%
%       col 1 (narrow):  baseline, zoomed 10 s window
%       col 2 (wide):    baseline, entire recording
%       col 3 (narrow):  treatment, zoomed 10 s window
%       col 4 (wide):    treatment, entire recording
%
%   Column widths are ~1:4.5 (zoom : full) following the
%   Chiappalone 2006 / Frega 2012 layout convention; the two sub-columns
%   of each condition sit directly adjacent so the eye reads them as a
%   single panel. The group titles "Baseline" and "<Treatment>" span
%   both sub-columns of their condition.
%
%   Each row is a single channel. The y-axis is scaled per-row to the
%   larger of (baseline, treatment) 99.5 % of |x|, rounded to a nice µV
%   value. A vertical scale bar with the µV value sits at the far right
%   of each row. Time scale bars (10 s, 60 s) are drawn under the
%   bottom row under their respective columns.
%
%   This is Panel A of the new Figure 2 (DOI) / Figure 3 (Ketanserin)
%   composites; Panel B is produced by fig_rate_panel_raster.m.
%
%   Options (Name-Value):
%     'pairIndex'     which pair from get_pairs_and_labels (default: mid)
%     'channels'      explicit (1 x N) channel list (default: auto)
%     'nChannels'     number of rows when channels is auto (default 8)
%     'zoomSec'       zoom window width in seconds (default 10)
%     'zoomStartSec' zoom start time in seconds (default 30)
%     'decimFactor'   decimation for the "entire recording" columns
%                     (default 4; preserves spike-band detail at 77 mm width)
%     'outName'       override the output filename stem (default
%                     '<study>_rate_panel_A_traces')

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pairIndex',    []);
    addParameter(p, 'channels',     []);
    addParameter(p, 'nChannels',    8);
    addParameter(p, 'zoomSec',      10);
    addParameter(p, 'zoomStartSec', 30);
    addParameter(p, 'decimFactor',  4);
    addParameter(p, 'outName',      '');
    addParameter(p, 'outDir',       fullfile(output_path(cfg, study, 'rates', ''), 'traces'));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_rate_panel_traces:NoPairs', 'No pairs for study %s.', study);
    end
    if isempty(opt.pairIndex)
        pairIdx = max(1, round(numel(pairs) / 2));
    else
        pairIdx = min(opt.pairIndex, numel(pairs));
    end
    pair = pairs(pairIdx);

    fprintf('fig_rate_panel_traces(%s): pair %d of %d\n  baseline: %s\n  treatment: %s\n', ...
        study, pairIdx, numel(pairs), pair.baseline, pair.treatment);

    % --- Choose channels ------------------------------------------------
    bCache = load_cache(pair.baseline,  cfg);
    tCache = load_cache(pair.treatment, cfg);
    if isempty(opt.channels)
        chans = select_representative_channels( ...
            bCache.spikeRates, tCache.spikeRates, opt.nChannels);
    else
        chans = opt.channels(:)';
    end
    nCh = numel(chans);
    fprintf('  selected channels: %s\n', mat2str(chans));

    % --- TDT SDK on path so TDTbin2mat resolves ------------------------
    addpath(genpath(cfg.paths.tdt_sdk));

    % --- Load + filter selected channels for both recordings ----------
    bBlock = dataset_path(pair.baseline,  cfg);
    tBlock = dataset_path(pair.treatment, cfg);

    bTraces = cell(1, nCh); tTraces = cell(1, nCh);
    bFs     = NaN;          tFs     = NaN;
    bDur    = NaN;          tDur    = NaN;

    fprintf('  loading %d channels per recording (~3-5 s each)...\n', nCh);
    t0 = tic;
    for i = 1:nCh
        ch = chans(i);
        try
            [x, fs, dur, ~] = load_and_filter(bBlock, ch, cfg);
            bTraces{i} = single(x);
            if isnan(bFs);  bFs  = fs;  end
            if isnan(bDur); bDur = dur; end
        catch ME
            warning('fig_rate_panel_traces:BaselineLoad', ...
                'ch %d: %s', ch, ME.message);
            bTraces{i} = single(zeros(1, 1));
        end
        try
            [x, fs, dur, ~] = load_and_filter(tBlock, ch, cfg);
            tTraces{i} = single(x);
            if isnan(tFs);  tFs  = fs;  end
            if isnan(tDur); tDur = dur; end
        catch ME
            warning('fig_rate_panel_traces:TreatmentLoad', ...
                'ch %d: %s', ch, ME.message);
            tTraces{i} = single(zeros(1, 1));
        end
        fprintf('    ch %2d  (%d/%d)  [%.0f s]\n', ch, i, nCh, toc(t0));
    end

    % --- Per-row amplitude scale (joint baseline+treatment) -----------
    % For the trace-display amplitude we use the true max(|x|) across
    % both conditions so that no spike is clipped.  A noise-floor lower
    % bound (8 * robust SD) protects channels whose max is negligible.
    rowScaleUv = zeros(1, nCh);
    for i = 1:nCh
        bx = abs(double(bTraces{i})) * 1e6;
        tx = abs(double(tTraces{i})) * 1e6;
        peakB = max(bx);
        peakT = max(tx);
        % Noise floor estimate (MAD)
        nB = mad(double(bTraces{i}), 1) * 1e6 / 0.6745;
        nT = mad(double(tTraces{i}), 1) * 1e6 / 0.6745;
        raw = max([peakB, peakT, 8 * nB, 8 * nT]);
        if raw <= 0 || ~isfinite(raw); raw = 50; end
        rowScaleUv(i) = nice_round_uv(raw * 1.05);
    end

    % =========================================================================
    % Figure layout (manual axes positioning, normalized to [0,1])
    % =========================================================================
    % Widths (zoom : full = 1 : 4.5), plus margins and a condition gap.
    leftMargin   = 0.03;
    rightMargin  = 0.07;   % room for µV scale-bar labels
    topMargin    = 0.06;   % room for group titles
    bottomMargin = 0.10;   % room for time scale bars + labels
    condGap      = 0.04;   % gap between baseline and treatment groups
    zoomFullGap  = 0.005;  % tiny gap so the eye separates zoom from full
    rowGap       = 0.004;

    availW = 1 - leftMargin - rightMargin - condGap - 2*zoomFullGap;
    wUnit  = availW / (2 * (1 + 4.5));
    wZoom  = wUnit;
    wFull  = wUnit * 4.5;

    xB1 = leftMargin;
    xB2 = xB1 + wZoom + zoomFullGap;
    xT1 = xB2 + wFull + condGap;
    xT2 = xT1 + wZoom + zoomFullGap;

    availH = 1 - topMargin - bottomMargin;
    rowH   = (availH - (nCh - 1) * rowGap) / nCh;

    % Panel sized for the LEFT 42% slot in a 183 mm composite (77 mm wide).
    % At this slot width, text and line sizes survive without down-scaling.
    % Height reduced proportionally; individual spike waveforms remain
    % visible (Chiappalone-style).
    fig = create_panel_figure(7.7, 5.0);

    % --- Group titles ---------------------------------------------------
    hdrB = xB1 + (wZoom + zoomFullGap + wFull) / 2;
    hdrT = xT1 + (wZoom + zoomFullGap + wFull) / 2;
    annotation(fig, 'textbox', [hdrB-0.1, 0.955, 0.2, 0.035], ...
        'String', labels.baseline, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');
    annotation(fig, 'textbox', [hdrT-0.1, 0.955, 0.2, 0.035], ...
        'String', labels.treatment, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold');

    % --- Trace rows -----------------------------------------------------
    zoomStart = opt.zoomStartSec;
    zoomEnd   = zoomStart + opt.zoomSec;

    for i = 1:nCh
        y = 1 - topMargin - i * rowH - (i - 1) * rowGap;
        yLim = rowScaleUv(i) * [-1 1];

        ax = axes('Position', [xB1 y wZoom rowH]); %#ok<LAXES>
        plot_trace(ax, bTraces{i}, bFs, [zoomStart zoomEnd], 1,               yLim);

        ax = axes('Position', [xB2 y wFull rowH]); %#ok<LAXES>
        plot_trace(ax, bTraces{i}, bFs, [0 bDur],             opt.decimFactor, yLim);

        ax = axes('Position', [xT1 y wZoom rowH]); %#ok<LAXES>
        plot_trace(ax, tTraces{i}, tFs, [zoomStart zoomEnd], 1,               yLim);

        ax = axes('Position', [xT2 y wFull rowH]); %#ok<LAXES>
        plot_trace(ax, tTraces{i}, tFs, [0 tDur],             opt.decimFactor, yLim);

        % µV scale bar on the far right (vertical tick + label).
        scaleX = xT2 + wFull + 0.006;
        annotation(fig, 'line', [scaleX scaleX], [y + rowH*0.25, y + rowH*0.75], ...
            'Color', 'k', 'LineWidth', 0.6);
        annotation(fig, 'textbox', [scaleX + 0.006, y + rowH*0.15, 0.09, rowH*0.7], ...
            'String', sprintf('%g \\muV', rowScaleUv(i)), ...
            'EdgeColor', 'none', 'FontSize', 6, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end

    % --- Time scale bars underneath the bottom row ----------------------
    scaleY = 0.055;   % absolute y in figure coords (fits above textbox)

    % 10 s bars under the two zoom columns.
    barZoomFracX = opt.zoomSec / opt.zoomSec;  % full width of the zoom axes
    barZoomW = wZoom * barZoomFracX * 0.5;     % bar is half of the axes width

    annotation(fig, 'line', [xB1, xB1 + barZoomW], [scaleY scaleY], ...
        'Color', 'k', 'LineWidth', 0.6);
    annotation(fig, 'textbox', [xB1, 0.008, wZoom, 0.04], ...
        'String', '10 s', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'FontSize', 6, 'FontWeight', 'bold');

    annotation(fig, 'line', [xT1, xT1 + barZoomW], [scaleY scaleY], ...
        'Color', 'k', 'LineWidth', 0.6);
    annotation(fig, 'textbox', [xT1, 0.008, wZoom, 0.04], ...
        'String', '10 s', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'FontSize', 6, 'FontWeight', 'bold');

    % 60 s bars under the two full columns (fraction of full duration).
    if bDur > 0
        barFullW = wFull * (60 / bDur);
    else
        barFullW = wFull * 0.1;
    end

    annotation(fig, 'line', [xB2, xB2 + barFullW], [scaleY scaleY], ...
        'Color', 'k', 'LineWidth', 0.6);
    annotation(fig, 'textbox', [xB2, 0.008, wFull, 0.04], ...
        'String', '60 s', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'FontSize', 6, 'FontWeight', 'bold');

    annotation(fig, 'line', [xT2, xT2 + barFullW], [scaleY scaleY], ...
        'Color', 'k', 'LineWidth', 0.6);
    annotation(fig, 'textbox', [xT2, 0.008, wFull, 0.04], ...
        'String', '60 s', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'FontSize', 6, 'FontWeight', 'bold');

    % --- Nature/NPP styling -----------------------------------------------
    apply_nature_style(fig);

    % --- Save -----------------------------------------------------------
    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    if isempty(opt.outName)
        outBase = sprintf('%s_rate_panel_A_traces', study);
    else
        outBase = opt.outName;
    end
    outFile = fullfile(opt.outDir, [outBase '.png']);
    save_figure(fig, outFile);
    close(fig);

    fprintf('fig_rate_panel_traces(%s): saved %s\n', study, outBase);
end

% =========================================================================
function plot_trace(ax, traceVolts, fs, tRangeSec, decim, yLimUv)
    nSamp = numel(traceVolts);
    if nSamp < 2 || isnan(fs)
        axis(ax, 'off'); return;
    end
    iStart = max(1, round(tRangeSec(1) * fs) + 1);
    iEnd   = min(nSamp, round(tRangeSec(2) * fs));
    if iEnd <= iStart
        axis(ax, 'off'); return;
    end
    seg = traceVolts(iStart:iEnd);
    if decim > 1
        seg = seg(1:decim:end);
        segFs = fs / decim;
    else
        segFs = fs;
    end
    t = (0:numel(seg)-1) / segFs + tRangeSec(1);
    plot(ax, t, double(seg) * 1e6, 'k', 'LineWidth', 0.3);
    xlim(ax, tRangeSec);
    ylim(ax, yLimUv);
    axis(ax, 'off');
end

% =========================================================================
function v = nice_round_uv(raw)
    if raw < 10
        v = 2 * ceil(raw / 2);
    elseif raw < 100
        v = 20 * ceil(raw / 20);
    elseif raw < 500
        v = 20 * ceil(raw / 20);
    elseif raw < 1000
        v = 40 * ceil(raw / 40);
    else
        v = 100 * ceil(raw / 100);
    end
end
