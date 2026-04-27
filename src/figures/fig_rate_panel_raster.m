function fig_rate_panel_raster(study, varargin)
%FIG_RATE_PANEL_RASTER Full-recording spike-time raster for one pair.
%
%   FIG_RATE_PANEL_RASTER(study) draws a two-column raster plot for a
%   representative pair from STUDY ('doi' or 'ket'):
%
%       left column  - baseline,  all channels, entire recording
%       right column - treatment, all channels, entire recording
%
%   Y axis is channel index; X axis is time in seconds; spikes are
%   drawn as short vertical ticks. This is Panel B of the new Figure 2
%   (DOI) / Figure 3 (Ketanserin) composites; Panel A is produced by
%   fig_rate_panel_traces.m.
%
%   FIG_RATE_PANEL_RASTER(study, 'Name', value, ...) options:
%
%     'pairIndex'  which pair from get_pairs_and_labels (default: mid)
%     'channels'   explicit channel list (default: cfg.channels.default)
%     'tickWidth'  line width of each spike tick (default 0.4)
%     'outName'    override output filename stem
%
% INPUTS:
%   study  -  'doi' | 'ket'
%
% OUTPUTS:
%   PDF + PNG sidecars written to output/fig{2,4}/panels/.
%
% See also: FIG_RATE_PANEL_TRACES, LOAD_PAIR_SPIKES, LOAD_CACHE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pairIndex', []);
    addParameter(p, 'channels',  cfg.channels.default);
    addParameter(p, 'tickWidth', 0.4);
    addParameter(p, 'outName',   '');
    addParameter(p, 'outDir',    fullfile(output_path(cfg, study, 'rates', ''), 'raster'));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_rate_panel_raster:NoPairs', 'No pairs for study %s.', study);
    end
    if isempty(opt.pairIndex)
        pairIdx = max(1, round(numel(pairs) / 2));
    else
        pairIdx = min(opt.pairIndex, numel(pairs));
    end
    pair = pairs(pairIdx);

    fprintf('fig_rate_panel_raster(%s): pair %d of %d\n  baseline: %s\n  treatment: %s\n', ...
        study, pairIdx, numel(pairs), pair.baseline, pair.treatment);

    % --- Load spike times from cache (no re-processing) ----------------
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, opt.channels, cfg);

    % Panel sized for the LEFT 42% slot in a 183 mm composite (77 mm wide).
    % Height reduced proportionally so text/lines are at final print size.
    fig = create_panel_figure(7.7, 3.5);
    tl = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Raster: %s  vs  %s   (pair %d)', ...
        labels.baseline, labels.treatment, pairIdx), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    axB = nexttile(tl, 1);
    plot_raster_column(axB, bSpikes, bMeta.durationSec, ...
        labels.baseline, opt.tickWidth);

    axT = nexttile(tl, 2);
    plot_raster_column(axT, tSpikes, tMeta.durationSec, ...
        labels.treatment, opt.tickWidth);

    % --- Nature/NPP styling -----------------------------------------------
    apply_nature_style(fig);

    % --- Save -----------------------------------------------------------
    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    if isempty(opt.outName)
        outBase = sprintf('%s_rate_panel_B_raster', study);
    else
        outBase = opt.outName;
    end
    outFile = fullfile(opt.outDir, [outBase '.png']);
    save_figure(fig, outFile);
    close(fig);

    fprintf('fig_rate_panel_raster(%s): saved %s\n', study, outBase);
end

% =========================================================================
function plot_raster_column(ax, spikeTimes, durSec, ttl, tickW)
    nCh = numel(spikeTimes);
    hold(ax, 'on');

    % Build one big NaN-separated XY vector so the entire raster is a
    % single line() call (dramatically faster than nCh * nSpikesPerCh
    % individual plot calls).
    totalSpikes = sum(cellfun(@numel, spikeTimes));
    Xall = nan(3 * totalSpikes, 1);
    Yall = nan(3 * totalSpikes, 1);
    idx = 1;
    for c = 1:nCh
        ts = spikeTimes{c};
        n = numel(ts);
        if n == 0; continue; end
        k = idx:(idx + 3*n - 1);
        Xall(k(1:3:end)) = ts(:);
        Xall(k(2:3:end)) = ts(:);
        Xall(k(3:3:end)) = NaN;
        Yall(k(1:3:end)) = c - 0.4;
        Yall(k(2:3:end)) = c + 0.4;
        Yall(k(3:3:end)) = NaN;
        idx = idx + 3*n;
    end
    line(ax, Xall, Yall, 'Color', 'k', 'LineWidth', tickW);

    xlim(ax, [0 max(durSec, 1)]);
    ylim(ax, [0.5 nCh + 0.5]);
    set(ax, 'YDir', 'normal');
    xlabel(ax, 'Time (s)', 'FontWeight', 'bold');
    ylabel(ax, 'Channel',   'FontWeight', 'bold');
    title(ax, ttl, 'FontWeight', 'bold');
    box(ax, 'on');
    set(ax, 'FontSize', 6);
    hold(ax, 'off');
end
