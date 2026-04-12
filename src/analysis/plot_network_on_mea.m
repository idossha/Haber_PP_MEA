function ax = plot_network_on_mea(adjacency, varargin)
%PLOT_NETWORK_ON_MEA Draw a connectivity graph over the 60MEA 8x8 geometry.
%
%   ax = PLOT_NETWORK_ON_MEA(adjacency) draws nodes at the spatial
%   positions of the 60 electrodes of the Multi Channel Systems
%   60MEA100/10iR layout, and edges between nodes where the weight is
%   above a threshold. Node colour / size can encode a per-channel scalar
%   (e.g. spike rate or degree), and edges can be coloured by weight or
%   by a signed difference (e.g. baseline vs treatment).
%
%   PLOT_NETWORK_ON_MEA(adjacency, 'Name', value, ...) with any of:
%
%     'parent'         axes handle (default: gca)
%     'nodeMetric'     (1 x nCh) per-channel scalar for node colour/size
%                      (default: node degree from the binarised adjacency)
%     'nodeCmap'       colormap name or matrix for nodes (default 'parula')
%     'nodeCLim'       [lo hi] for node colour mapping (default data range)
%     'nodeSizeRange'  [min max] marker area pts^2 (default [30 400])
%     'edgeMode'       'weight' | 'delta' (default 'weight')
%     'edgeThresholdPct'  keep the top P percent of edges by |value|
%                      (default 20). Set to [] to use the explicit
%                      'edgeThreshold' instead.
%     'edgeThreshold'  explicit numeric threshold; overrides the percent
%     'edgeCmap'       colormap for 'weight' mode (default 'gray')
%     'edgeDivCmap'    diverging colormap for 'delta' (default blue-white-red)
%     'edgeCLim'       [lo hi] colour limits; default symmetric for delta,
%                      data range for weight
%     'edgeAlphaRange' [min max] alpha mapping as function of |value|
%                      (default [0.15 0.9])
%     'edgeWidthRange' [min max] line width (default [0.5 3])
%     'title'          axis title (default '')
%     'showCorners'    logical, draw grey rectangles at the 4 missing
%                      corner positions (default false)
%
% INPUTS:
%   adjacency  -  nCh x nCh matrix (symmetric, weighted). For 'delta'
%                 mode pass (treatment - baseline). NaNs and the diagonal
%                 are ignored.
%
% OUTPUTS:
%   ax  -  The axes handle used for plotting (matches 'parent' if given).
%
% Notes:
%   - Layout is read from src/utils/mea60_layout.m (single source of
%     truth). Channels outside 1..60 are ignored.
%   - To draw baseline and treatment side-by-side, call this function
%     twice with different axes handles (see fig_connectivity_exemplar.m).
%
% See also: CONNECTIVITY_XCORR, NETWORK_METRICS, TOPOGRAPHICAL_MAP,
%           FIG_CONNECTIVITY_EXEMPLAR.

    p = inputParser;
    addRequired(p,  'adjacency', @(x) isnumeric(x) && ismatrix(x) && size(x,1) == size(x,2));
    addParameter(p, 'parent',           []);
    addParameter(p, 'nodeMetric',       []);
    addParameter(p, 'nodeCmap',         'parula');
    addParameter(p, 'nodeCLim',         []);
    addParameter(p, 'nodeSizeRange',    [30 400]);
    addParameter(p, 'edgeMode',         'weight', @(s) any(strcmpi(s, {'weight','delta'})));
    addParameter(p, 'edgeThresholdPct', 20, @(x) isempty(x) || (isscalar(x) && x >= 0 && x <= 100));
    addParameter(p, 'edgeThreshold',    [], @(x) isempty(x) || isscalar(x));
    addParameter(p, 'edgeCmap',         'gray');
    addParameter(p, 'edgeDivCmap',      []);
    addParameter(p, 'edgeCLim',         []);
    addParameter(p, 'edgeAlphaRange',   [0.15 0.9]);
    addParameter(p, 'edgeWidthRange',   [0.5 3]);
    addParameter(p, 'title',            '');
    addParameter(p, 'showCorners',      false, @(x) islogical(x) && isscalar(x));
    parse(p, adjacency, varargin{:});
    opt = p.Results;

    if isempty(opt.parent)
        ax = gca;
    else
        ax = opt.parent;
    end
    hold(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'none', 'YColor', 'none');
    axis(ax, 'equal');
    set(ax, 'XLim', [0.3 8.7], 'YLim', [0.3 8.7]);
    if ~isempty(opt.title)
        title(ax, opt.title);
    end

    layout = mea60_layout();
    nCh    = layout.nCh;

    % --- Grid frame and missing-corner markers (optional) ---
    if opt.showCorners
        for r = [1 8]
            for c = [1 8]
                rectangle(ax, ...
                    'Position', [c-0.3, (9-r)-0.3, 0.6, 0.6], ...
                    'EdgeColor', [0.85 0.85 0.85], 'LineStyle', ':');
            end
        end
    end

    % --- Node colour/size metric ---
    if isempty(opt.nodeMetric)
        % Default: node degree from the above-threshold adjacency
        A = adjacency;
        A(isnan(A)) = 0;
        A(1:nCh+1:end) = 0;
        tmpThreshold = edge_threshold_value(A, opt);
        bin = (A > tmpThreshold) & ~isnan(A);
        nodeMetric = sum(bin, 2);
    else
        nodeMetric = opt.nodeMetric(:);
        if numel(nodeMetric) < nCh
            nodeMetric(end+1:nCh) = NaN;
        end
    end

    % --- Edge selection and plotting ---
    A = adjacency;
    A(1:nCh+1:end) = NaN;
    upperMask = triu(true(nCh), 1);

    upperVals = A(upperMask);
    upperVals(isnan(upperVals)) = 0;

    thresholdValue = edge_threshold_value(A, opt);

    switch lower(opt.edgeMode)
        case 'weight'
            keepVals = upperVals;
            edgeKeep = keepVals >= thresholdValue;
        case 'delta'
            keepVals = upperVals;                % signed values
            edgeKeep = abs(keepVals) >= thresholdValue;
    end

    % Determine edge colour limits.
    if isempty(opt.edgeCLim)
        switch lower(opt.edgeMode)
            case 'weight'
                if any(edgeKeep)
                    edgeCLim = [min(keepVals(edgeKeep)), max(keepVals(edgeKeep))];
                else
                    edgeCLim = [0 1];
                end
            case 'delta'
                if any(edgeKeep)
                    maxAbs = max(abs(keepVals(edgeKeep)));
                    edgeCLim = [-maxAbs, maxAbs];
                else
                    edgeCLim = [-1 1];
                end
        end
    else
        edgeCLim = opt.edgeCLim;
    end

    % Build colormap for edges.
    switch lower(opt.edgeMode)
        case 'weight'
            edgeCmap = get_colormap(opt.edgeCmap);
        case 'delta'
            if isempty(opt.edgeDivCmap)
                edgeCmap = diverging_cmap(256);
            else
                edgeCmap = get_colormap(opt.edgeDivCmap);
            end
    end

    % Draw the edges in ascending order of |value| so strong ones land on top.
    [~, order] = sort(abs(keepVals));
    [iIdx, jIdx] = find(upperMask);
    iIdx = iIdx(order);
    jIdx = jIdx(order);
    keepValsSorted = keepVals(order);
    edgeKeepSorted = edgeKeep(order);

    for e = 1:numel(keepValsSorted)
        if ~edgeKeepSorted(e); continue; end
        v = keepValsSorted(e);
        ci = iIdx(e);
        cj = jIdx(e);
        x1 = layout.xCoord(ci);
        y1 = layout.yCoord(ci);
        x2 = layout.xCoord(cj);
        y2 = layout.yCoord(cj);

        t = (v - edgeCLim(1)) / max(edgeCLim(2) - edgeCLim(1), eps);
        t = min(max(t, 0), 1);
        rowIdx = 1 + round(t * (size(edgeCmap, 1) - 1));
        color = edgeCmap(rowIdx, :);

        aT = min(max(abs(v) / max(abs(edgeCLim)), 0), 1);
        alpha = opt.edgeAlphaRange(1) + ...
                (opt.edgeAlphaRange(2) - opt.edgeAlphaRange(1)) * aT;
        width = opt.edgeWidthRange(1) + ...
                (opt.edgeWidthRange(2) - opt.edgeWidthRange(1)) * aT;

        line(ax, [x1 x2], [y1 y2], ...
             'Color', [color alpha], ...
             'LineWidth', max(width, 0.01));
    end

    % --- Nodes ---
    nodeCmap = get_colormap(opt.nodeCmap);
    finiteNode = nodeMetric(isfinite(nodeMetric));
    if isempty(opt.nodeCLim)
        if isempty(finiteNode)
            nodeCLim = [0 1];
        else
            nodeCLim = [min(finiteNode), max(finiteNode)];
            if diff(nodeCLim) == 0
                nodeCLim = nodeCLim + [-1 1];
            end
        end
    else
        nodeCLim = opt.nodeCLim;
    end

    for c = 1:nCh
        x = layout.xCoord(c);
        y = layout.yCoord(c);
        v = nodeMetric(c);
        if isnan(v)
            faceColor = [0.8 0.8 0.8];
            sz = opt.nodeSizeRange(1);
        else
            tNode = (v - nodeCLim(1)) / max(nodeCLim(2) - nodeCLim(1), eps);
            tNode = min(max(tNode, 0), 1);
            rowIdx = 1 + round(tNode * (size(nodeCmap, 1) - 1));
            faceColor = nodeCmap(rowIdx, :);
            sz = opt.nodeSizeRange(1) + ...
                 (opt.nodeSizeRange(2) - opt.nodeSizeRange(1)) * tNode;
        end
        scatter(ax, x, y, sz, faceColor, 'filled', ...
                'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 0.75);
    end

    hold(ax, 'off');
end

% =========================================================================
function t = edge_threshold_value(A, opt)
% Resolve the explicit 'edgeThreshold' if given; otherwise use the
% 'edgeThresholdPct' percentile of the |A| upper triangle.
    if ~isempty(opt.edgeThreshold)
        t = opt.edgeThreshold;
        return;
    end
    n = size(A, 1);
    upperMask = triu(true(n), 1);
    vals = A(upperMask);
    vals = vals(~isnan(vals));
    if strcmpi(opt.edgeMode, 'delta')
        vals = abs(vals);
    end
    if isempty(vals)
        t = 0;
        return;
    end
    pct = 100 - opt.edgeThresholdPct;
    t = prctile(vals, pct);
end

% =========================================================================
function cmap = get_colormap(name)
% Resolve a colormap name or accept an explicit Nx3 matrix.
    if ischar(name) || isstring(name)
        fcn = str2func(char(name));
        cmap = fcn(256);
    elseif isnumeric(name) && size(name, 2) == 3
        cmap = name;
    else
        error('plot_network_on_mea:badCmap', ...
            'Colormap must be a name or an Nx3 matrix.');
    end
end

% =========================================================================
function cmap = diverging_cmap(n)
% Simple blue-white-red diverging colormap for delta panels.
    half = floor(n / 2);
    top  = n - half;
    blue = [linspace(0.11,1,half).', linspace(0.30,1,half).', linspace(0.60,1,half).'];
    red  = [linspace(1,0.75,top).',  linspace(1,0.10,top).',  linspace(1,0.15,top).'];
    cmap = [blue; red];
end
