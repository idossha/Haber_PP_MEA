function fig_rate_change_bar(metric, study, varargin)
%FIG_RATE_CHANGE_BAR Per-recording silenced/decreased/increased boxplots.
%
%   FIG_RATE_CHANGE_BAR(metric, study) categorises baseline-active channels
%   into Silenced (treatment == 0), Decreased (still active, lower), or
%   Increased, then renders a three-box boxplot of per-recording percentages
%   and prints a chi-square goodness-of-fit (50/50) on the increased vs
%   decreased counts (silenced excluded). Replaces both
%   figures_spike_rate_change_bar.m and figures_burst_rate_change_bar.m.
%
% INPUTS:
%   metric  -  'spike' or 'burst' (which cache field to categorise).
%   study   -  'doi' or 'ket'.
%
% Name-value options:
%   'channels'              -  default cfg.channels.default
%   'rateEqualRelTol'       -  default 1e-6 (unchanged tolerance)
%   'ignoreSilentChannels'  -  burst metric only, default true
%   'minRateThreshold'      -  burst metric only, default 0 (bursts/min)
%
% OUTPUTS:
%   PNG written to output/fig{2,4}/panels/<study>_<metric>_rate_change_categories.png
%   plus chi-square text on stdout.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'metric', @(s) any(strcmpi(s, {'spike','burst'})));
    addRequired(p,  'study',  @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',              cfg.channels.default);
    addParameter(p, 'rateEqualRelTol',       1e-6);
    addParameter(p, 'ignoreSilentChannels',  true);
    addParameter(p, 'minRateThreshold',      0);
    parse(p, metric, study, varargin{:});
    opt    = p.Results;
    metric = lower(opt.metric);
    study  = lower(opt.study);

    pairs = get_pairs_and_labels(cfg, study);

    cacheField = 'spikeRates';
    if strcmp(metric, 'burst')
        cacheField = 'burstRates';
    end

    % Per-recording category counts (totals are accumulated for chi-square).
    nSilenced  = 0;
    nDecreased = 0;
    nIncreased = 0;
    nUnchanged = 0;
    pctSil = [];
    pctDec = [];
    pctInc = [];

    for k = 1:numel(pairs)
        [loadedB, loadedT] = load_pair_cache(pairs(k), cfg);
        ratesB = loadedB.(cacheField)(:);
        ratesT = loadedT.(cacheField)(:);
        chB    = loadedB.channelsUsed(:)';
        chT    = loadedT.channelsUsed(:)';

        [~, idxB] = ismember(opt.channels, chB);
        [~, idxT] = ismember(opt.channels, chT);
        valid = (idxB > 0) & (idxT > 0);
        if ~any(valid); continue; end

        b = ratesB(idxB(valid));
        t = ratesT(idxT(valid));
        ok = ~isnan(b) & ~isnan(t);
        b = b(ok); t = t(ok);

        if strcmp(metric, 'burst') && opt.ignoreSilentChannels
            if opt.minRateThreshold == 0
                nonSilent = ~(b == 0 & t == 0);
            else
                nonSilent = ~(b < opt.minRateThreshold & t < opt.minRateThreshold);
            end
            b = b(nonSilent);
            t = t(nonSilent);
        end
        if isempty(b); continue; end

        % Baseline-active only.
        use = b > 0;
        b2 = b(use); t2 = t(use);
        if isempty(b2); continue; end

        silent = (t2 == 0);
        nS = sum(silent);
        nSilenced = nSilenced + nS;

        active = ~silent;
        b3 = b2(active);
        t3 = t2(active);
        if isempty(b3)
            nD = 0; nI = 0; nU = 0;
        else
            d = t3 - b3;
            thresh = opt.rateEqualRelTol .* max(b3, eps);
            unch = abs(d) <= thresh;
            inc  = d > 0 & ~unch;
            dec  = d < 0 & ~unch;
            nD = sum(dec);
            nI = sum(inc);
            nU = sum(unch);
        end
        nIncreased = nIncreased + nI;
        nDecreased = nDecreased + nD;
        nUnchanged = nUnchanged + nU;

        nDenom = nS + nD + nI;
        if nDenom == 0; continue; end
        pctSil(end+1) = 100 * nS / nDenom; %#ok<AGROW>
        pctDec(end+1) = 100 * nD / nDenom; %#ok<AGROW>
        pctInc(end+1) = 100 * nI / nDenom; %#ok<AGROW>
    end

    nClassified = nSilenced + nDecreased + nIncreased + nUnchanged;
    if nClassified == 0
        error('fig_rate_change_bar:NoData', ...
            'No baseline-active channels found (need %s rate > 0).', metric);
    end

    fprintf('fig_rate_change_bar(%s, %s): n=%d  (silenced %d, decreased %d, increased %d, unchanged %d)\n', ...
        metric, study, nClassified, nSilenced, nDecreased, nIncreased, nUnchanged);

    chi2Result = compute_increase_decrease_chi2(nIncreased, nDecreased);
    print_increase_decrease_chi2(chi2Result, nIncreased, nDecreased);

    if isempty(pctSil)
        error('fig_rate_change_bar:NoData', 'No recordings with non-empty denominator.');
    end

    panelDir = output_path(cfg, study, 'rates', 'panels');
    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end
    outFile = fullfile(panelDir, ...
        sprintf('%s_%s_rate_change_categories.png', study, metric));

    [~, labels] = get_pairs_and_labels(cfg, study);
    if strcmp(metric, 'spike')
        metricName = 'Firing rate';
    else
        metricName = 'Burst rate';
    end
    titleText = '';

    plot_change_category_box(pctSil, pctDec, pctInc, nUnchanged, outFile, ...
        'title',  titleText, ...
        'yLabel', '% of baseline-active channels', ...
        'chi2',   chi2Result);
    fprintf('  saved: %s\n', outFile);

    % --- Numeric sidecar ------------------------------------------------
    stats = struct( ...
        'study',             study, ...
        'metric',            metric, ...
        'n_recordings',      numel(pctSil), ...
        'n_silenced',        nSilenced, ...
        'n_decreased',       nDecreased, ...
        'n_increased',       nIncreased, ...
        'n_unchanged',       nUnchanged, ...
        'median_pct_silenced',  median(pctSil), ...
        'median_pct_decreased', median(pctDec), ...
        'median_pct_increased', median(pctInc), ...
        'mean_pct_silenced',    mean(pctSil), ...
        'mean_pct_decreased',   mean(pctDec), ...
        'mean_pct_increased',   mean(pctInc), ...
        'chi2_statistic',    chi2Result.stat, ...
        'chi2_df',           chi2Result.df, ...
        'chi2_p',            chi2Result.p, ...
        'pct_increase_of_directional',  chi2Result.pctInc, ...
        'pct_decrease_of_directional',  chi2Result.pctDec, ...
        'figure_file',       outFile);
    export_figure_stats(stats, fullfile(statsDir, ...
        sprintf('%s_%s_rate_change_stats', study, metric)));
end

% =========================================================================
function plot_change_category_box(pctSil, pctDec, pctInc, nUnchanged, outFile, varargin)
    p = inputParser;
    addRequired(p,  'pctSil');
    addRequired(p,  'pctDec');
    addRequired(p,  'pctInc');
    addRequired(p,  'nUnchanged');
    addRequired(p,  'outFile');
    addParameter(p, 'title',  '');
    addParameter(p, 'yLabel', '% of baseline-active channels');
    addParameter(p, 'chi2',   struct());
    parse(p, pctSil, pctDec, pctInc, nUnchanged, outFile, varargin{:});
    opt = p.Results;

    colorSilenced = [0.38, 0.38, 0.38];
    colorDecrease = [0.82, 0.58, 0.58];
    colorIncrease = [0.50, 0.72, 0.55];

    fig = create_panel_figure(5.8, 4.4);
    hold on;

    n = numel(pctSil);
    y = [pctSil(:); pctDec(:); pctInc(:)];
    g = [ones(n,1); 2*ones(n,1); 3*ones(n,1)];
    boxplot(y, g, 'Positions', [1 2 3], ...
        'Labels', {'Silenced', 'Decreased', 'Increased'}, ...
        'Widths', 0.55, 'Symbol', '');

    hBox = findobj(gca, 'Tag', 'Box');
    nx = zeros(numel(hBox), 1);
    for jb = 1:numel(hBox)
        xd = get(hBox(jb), 'XData');
        nx(jb) = mean(xd(:));
    end
    [~, sortIdx] = sort(nx);
    hBoxSorted = hBox(sortIdx);
    faceColors    = [colorSilenced; colorDecrease; colorIncrease];
    outlineColors = [0.22 0.22 0.22; 0.55 0.35 0.35; 0.28 0.45 0.28];
    for jb = 1:numel(hBoxSorted)
        xd = get(hBoxSorted(jb), 'XData');
        yd = get(hBoxSorted(jb), 'YData');
        hp = patch(xd, yd, faceColors(jb,:), 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        uistack(hp, 'bottom');
        set(hBoxSorted(jb), 'Color', outlineColors(jb,:), 'LineWidth', 1.2);
    end
    hMed = findobj(gca, 'Tag', 'Median');
    set(hMed, 'Color', [0 0 0], 'LineWidth', 1.5);
    hWh = [findobj(gca, 'Tag', 'Upper Whisker'); findobj(gca, 'Tag', 'Lower Whisker')];
    if ~isempty(hWh)
        set(hWh, 'Color', [0.3 0.3 0.3], 'LineStyle', '-', 'LineWidth', 1.2);
    end
    hCap = [findobj(gca, 'Tag', 'Upper Adjacent Value'); findobj(gca, 'Tag', 'Lower Adjacent Value')];
    for k = 1:numel(hCap)
        xd = get(hCap(k), 'XData');
        xm = (xd(1) + xd(2)) / 2;
        d  = (xd(2) - xd(1)) * 0.5;
        set(hCap(k), 'XData', [xm - d/2, xm + d/2], 'LineWidth', 1.2);
    end

    meanVals = [mean(pctSil), mean(pctDec), mean(pctInc)];
    for jb = 1:3
        xd = get(hBoxSorted(jb), 'XData');
        plot([min(xd(:)), max(xd(:))], [meanVals(jb), meanVals(jb)], ...
            '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
    end

    box off;

    ylabel(opt.yLabel);
    if ~isempty(opt.title)
        title(opt.title, 'FontSize', 15, 'FontName', 'Arial', 'FontWeight', 'bold');
    end

    xlim([0.5, 3.5]);
    ylim([0, 100 * 1.05]);
    yTickStep = nice_tick_step(100 * 1.05, 8);
    yticks(0:yTickStep:ceil(100 * 1.05 / yTickStep) * yTickStep);

    % Chi-square result annotation (top-right).
    if isfield(opt.chi2, 'p') && ~isempty(opt.chi2) && isfield(opt.chi2, 'stat')
        chi2Str = sprintf('\\chi^2(%d) = %.2f,  p = %s', ...
            opt.chi2.df, opt.chi2.stat, format_p_value(opt.chi2.p));
        annotation(fig, 'textbox', [0.18, 0.89, 0.80, 0.05], ...
            'String', chi2Str, ...
            'FontSize', 14, 'FontName', 'Arial', 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    end

    if nUnchanged > 0
        annotation(fig, 'textbox', [0.18, 0.02, 0.80, 0.03], ...
            'String', sprintf('Unchanged channels (within tolerance): %d', nUnchanged), ...
            'FontSize', 11, 'FontName', 'Arial', 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end

    % --- Nature/NPP styling -----------------------------------------------
    apply_nature_style(fig);

    save_figure(fig, outFile);
    close(fig);
end

% -------------------------------------------------------------------------
function s = format_p_value(p)
    if isnan(p); s = 'n/a'; return; end
    if p < 0.001
        s = '<0.001';
    else
        s = sprintf('%.3f', p);
    end
end

% =========================================================================
function res = compute_increase_decrease_chi2(nIncreased, nDecreased)
% Goodness-of-fit of (increased, decreased) counts against a 50/50 null.
    res = struct('stat', NaN, 'df', 1, 'p', NaN, ...
                 'nInc', nIncreased, 'nDec', nDecreased, ...
                 'pctInc', NaN, 'pctDec', NaN);
    N = nIncreased + nDecreased;
    if N == 0; return; end
    obs = [nIncreased, nDecreased];
    expCounts = [N/2, N/2];
    res.stat = sum((obs - expCounts).^2 ./ expCounts);
    try
        res.p = 1 - chi2cdf(res.stat, res.df);
    catch
        res.p = NaN;
    end
    res.pctInc = 100 * nIncreased / N;
    res.pctDec = 100 * nDecreased / N;
end

% -------------------------------------------------------------------------
function print_increase_decrease_chi2(res, nIncreased, nDecreased)
    if isnan(res.stat)
        fprintf('  Chi-square: no directional change observations.\n');
        return;
    end
    if res.p < 0.001
        pStr = '<0.001';
    elseif res.p < 0.01
        pStr = sprintf('%.3f', res.p);
    else
        pStr = sprintf('%.3g', res.p);
    end
    fprintf('  Increase vs Decrease (excluding silenced): %.1f%% (%d) vs %.1f%% (%d)\n', ...
        res.pctInc, nIncreased, res.pctDec, nDecreased);
    fprintf('  Chi-square goodness-of-fit vs 50/50: chi2(%d) = %.4g, p = %s\n', ...
        res.df, res.stat, pStr);
    if ~isnan(res.p) && res.p < 0.05
        fprintf('  Deviation from 50/50 is significant (p < 0.05).\n');
    else
        fprintf('  Deviation from 50/50 is not significant (p >= 0.05).\n');
    end
end
