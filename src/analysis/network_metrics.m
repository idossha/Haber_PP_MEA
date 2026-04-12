function metrics = network_metrics(adjacency, varargin)
%NETWORK_METRICS Graph-level summaries from a connectivity adjacency matrix.
%
%   metrics = NETWORK_METRICS(adjacency) treats ADJACENCY as a symmetric,
%   undirected weighted graph (NaNs and the diagonal are ignored), and
%   returns a struct with the following fields:
%
%     Basic weight and degree statistics:
%       .meanWeight         -  mean of upper-triangular off-diagonal entries
%       .stdWeight          -  std of the same
%       .density            -  fraction of possible edges above threshold
%       .threshold          -  binarisation threshold actually used
%       .nNodes             -  number of nodes
%       .nEdges             -  number of edges above threshold
%       .degreeMean         -  mean node degree (binarised)
%       .degreeStd          -  std node degree (binarised)
%
%     Graph-theoretic summaries (binarised at threshold):
%       .clusteringMean     -  mean local clustering coefficient
%                              (Watts & Strogatz 1998)
%       .meanShortestPath   -  mean shortest-path length on the largest
%                              connected component (BFS)
%       .globalEfficiency   -  mean of 1/d(i,j) over all pairs (Latora &
%                              Marchiori 2001). 0 if disconnected pairs.
%       .modularity         -  Newman 2006 modularity Q for the partition
%                              returned by a greedy agglomerative scheme.
%                              NaN if the graph has < 2 edges.
%       .nCommunities       -  number of communities in that partition.
%       .smallWorldnessSigma -  sigma = (C/Crand) / (L/Lrand), using the
%                              mean of 'nullN' Erdos-Renyi random graphs
%                              with matched n and edge density
%                              (Humphries & Gurney 2008). NaN if graph
%                              is too small / disconnected for L.
%
%   NETWORK_METRICS(adjacency, 'threshold', t) sets the binarisation
%   threshold; default is the median of the off-diagonal upper triangle
%   ("threshold by data").
%
%   NETWORK_METRICS(adjacency, 'nullN', n) sets the number of random
%   null graphs for small-worldness (default 50, set 0 to skip).
%
%   NETWORK_METRICS(adjacency, 'nullSeed', s) fixes the RNG seed for the
%   null ensemble (default 1).
%
% INPUTS:
%   adjacency  -  Square symmetric matrix from connectivity_xcorr().
%
% Name-value options:
%   'threshold'  -  Scalar threshold; default = median(offDiagUpper).
%   'nullN'      -  Null ensemble size for small-worldness; default 50.
%   'nullSeed'   -  RNG seed for null ensemble; default 1.
%
% OUTPUTS:
%   metrics  -  Struct with fields above.
%
% References (cited in docs/METRICS.md):
%   - Rubinov & Sporns 2010, NeuroImage 52:1059-1069.
%   - Watts & Strogatz 1998, Nature 393:440-442.
%   - Latora & Marchiori 2001, Phys Rev Lett 87:198701 (global efficiency).
%   - Newman 2006, PNAS 103:8577-8582 (modularity).
%   - Humphries & Gurney 2008, PLoS ONE 3:e0002051 (small-worldness sigma).
%
% Notes:
%   - All metrics are computed on the binarised graph. Weighted equivalents
%     are listed under "Deviations from literature" in docs/METRICS.md.
%   - The modularity solver is a simple greedy merge (not Louvain). For
%     n = 60 it converges in well under a second and gives a partition
%     close to the Louvain optimum; swap in a Louvain implementation if
%     you ever need the best-possible Q.

    if size(adjacency, 1) ~= size(adjacency, 2)
        error('network_metrics:Args', 'adjacency must be square.');
    end

    p = inputParser;
    addParameter(p, 'threshold', [],   @(x) isempty(x) || isscalar(x));
    addParameter(p, 'nullN',     50,   @(x) isscalar(x) && x >= 0);
    addParameter(p, 'nullSeed',  1,    @(x) isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    nNodes = size(adjacency, 1);
    A = adjacency;
    A(1:nNodes+1:end) = NaN;        % drop diagonal

    upperMask = triu(true(nNodes), 1);
    upperVals = A(upperMask);
    upperVals = upperVals(~isnan(upperVals));

    metrics.meanWeight = mean(upperVals);
    metrics.stdWeight  = std(upperVals);

    if isempty(opt.threshold)
        threshold = median(upperVals);
    else
        threshold = opt.threshold;
    end
    metrics.threshold = threshold;

    binA = A > threshold;
    binA(isnan(A)) = false;
    binA = binA | binA';            % enforce symmetry
    binA(1:nNodes+1:end) = false;   % no self-loops

    nEdges    = sum(binA(upperMask));
    nPossible = sum(upperMask(:));
    metrics.density = nEdges / max(nPossible, 1);
    metrics.nEdges  = nEdges;
    metrics.nNodes  = nNodes;

    degree = sum(binA, 2);
    metrics.degreeMean = mean(degree);
    metrics.degreeStd  = std(degree);

    metrics.clusteringMean    = local_clustering(binA);
    metrics.meanShortestPath  = mean_shortest_path(binA);
    metrics.globalEfficiency  = global_efficiency(binA);

    [Q, comms] = greedy_modularity(binA);
    metrics.modularity   = Q;
    metrics.nCommunities = numel(unique(comms));

    % Small-worldness sigma vs random-null ensemble.
    if opt.nullN > 0 && metrics.density > 0 && ~isnan(metrics.clusteringMean)
        metrics.smallWorldnessSigma = small_worldness_sigma( ...
            binA, metrics.clusteringMean, metrics.meanShortestPath, ...
            opt.nullN, opt.nullSeed);
    else
        metrics.smallWorldnessSigma = NaN;
    end
end

% =========================================================================
function meanC = local_clustering(binA)
% Watts-Strogatz local clustering coefficient (binarised undirected).
    n = size(binA, 1);
    C = zeros(n, 1);
    for i = 1:n
        nbrs = find(binA(i, :));
        k = numel(nbrs);
        if k < 2
            continue;
        end
        sub = binA(nbrs, nbrs);
        triangles = sum(sub(:)) / 2;     % undirected
        C(i) = (2 * triangles) / (k * (k - 1));
    end
    meanC = mean(C);
end

% =========================================================================
function L = mean_shortest_path(binA)
% Mean shortest-path length on the largest connected component (BFS).
    n = size(binA, 1);
    if n == 0
        L = NaN;
        return;
    end

    comps = connected_components(binA);
    compSizes = histcounts(comps, 1:max(comps)+1);
    [~, largest] = max(compSizes);
    nodes = find(comps == largest);
    m = numel(nodes);
    if m < 2
        L = NaN;
        return;
    end

    sub = binA(nodes, nodes);
    distSum = 0;
    pairCount = 0;
    for i = 1:m
        d = bfs_distances(sub, i);
        finiteD = d(~isinf(d));
        finiteD = finiteD(finiteD > 0);
        distSum   = distSum + sum(finiteD);
        pairCount = pairCount + numel(finiteD);
    end
    if pairCount == 0
        L = NaN;
    else
        L = distSum / pairCount;
    end
end

% =========================================================================
function E = global_efficiency(binA)
% Latora & Marchiori 2001 global efficiency = mean of 1 / d(i,j) over all
% ordered pairs (i != j). Disconnected pairs contribute 0.
    n = size(binA, 1);
    if n < 2
        E = NaN;
        return;
    end
    invDistSum = 0;
    nPairs = 0;
    for i = 1:n
        d = bfs_distances(binA, i);
        for j = 1:n
            if j == i; continue; end
            nPairs = nPairs + 1;
            if isfinite(d(j)) && d(j) > 0
                invDistSum = invDistSum + 1 / d(j);
            end
        end
    end
    if nPairs == 0
        E = NaN;
    else
        E = invDistSum / nPairs;
    end
end

% =========================================================================
function [Q, comm] = greedy_modularity(binA)
% Greedy agglomerative modularity (merges the pair of communities whose
% merge yields the largest Delta Q until no merge improves Q). Start with
% every node in its own community. Good enough for n ~ 60.
    n = size(binA, 1);
    m = sum(binA(:)) / 2;
    if m <= 0
        Q = NaN;
        comm = 1:n;
        return;
    end

    k = sum(binA, 2);
    comm = (1:n)';          % node -> community id
    commActive = true(n, 1);

    % Precompute connectivity between communities.
    % For n = 60 we can just update Q on each merge by brute force.
    Q = modularity_Q(binA, comm, k, m);
    improved = true;

    while improved
        improved = false;
        bestDelta = 0;
        bestPair  = [];
        activeIds = find(commActive);
        for a = 1:numel(activeIds)
            ca = activeIds(a);
            for b = a+1:numel(activeIds)
                cb = activeIds(b);
                testComm = comm;
                testComm(testComm == cb) = ca;
                Qtest = modularity_Q(binA, testComm, k, m);
                delta = Qtest - Q;
                if delta > bestDelta
                    bestDelta = delta;
                    bestPair  = [ca, cb];
                end
            end
        end
        if ~isempty(bestPair)
            comm(comm == bestPair(2)) = bestPair(1);
            commActive(bestPair(2)) = false;
            Q = Q + bestDelta;
            improved = true;
        end
    end

    % Relabel communities 1..nComm for cleanliness.
    [~, ~, comm] = unique(comm);
end

% -------------------------------------------------------------------------
function Q = modularity_Q(binA, comm, k, m)
% Newman 2006 Q = (1 / (2m)) * sum_ij [A_ij - k_i k_j / (2m)] * delta(c_i, c_j)
    n = size(binA, 1);
    sameComm = bsxfun(@eq, comm(:), comm(:).');
    expected = (k * k.') / (2 * m);
    Q = sum(sum((binA - expected) .* sameComm)) / (2 * m);
end

% =========================================================================
function sigma = small_worldness_sigma(binA, Cobs, Lobs, nullN, nullSeed)
% sigma = (Cobs / Crand) / (Lobs / Lrand)
% Crand and Lrand are the mean over an Erdos-Renyi ensemble matched on
% n and m.
    if isnan(Lobs)
        sigma = NaN;
        return;
    end
    n = size(binA, 1);
    m = sum(binA(:)) / 2;
    if m <= 0
        sigma = NaN;
        return;
    end

    rng(nullSeed, 'twister');
    Cs = nan(nullN, 1);
    Ls = nan(nullN, 1);
    upperMask = triu(true(n), 1);
    upperIdx  = find(upperMask);
    nPossible = numel(upperIdx);
    for r = 1:nullN
        sel = randperm(nPossible, m);
        B   = false(n);
        B(upperIdx(sel)) = true;
        B = B | B';
        Cs(r) = local_clustering(B);
        Ls(r) = mean_shortest_path(B);
    end
    Crand = mean(Cs(~isnan(Cs)));
    Lrand = mean(Ls(~isnan(Ls)));

    if isnan(Crand) || Crand == 0 || isnan(Lrand) || Lrand == 0
        sigma = NaN;
        return;
    end
    sigma = (Cobs / Crand) / (Lobs / Lrand);
end

% =========================================================================
function comps = connected_components(binA)
% Flood-fill connected components on the binarised adjacency.
    n = size(binA, 1);
    comps = zeros(n, 1);
    nextId = 0;
    for s = 1:n
        if comps(s) > 0; continue; end
        nextId = nextId + 1;
        stack = s;
        while ~isempty(stack)
            v = stack(end); stack(end) = [];
            if comps(v) > 0; continue; end
            comps(v) = nextId;
            nbrs = find(binA(v, :));
            stack = [stack, nbrs(comps(nbrs) == 0)]; %#ok<AGROW>
        end
    end
end

% =========================================================================
function d = bfs_distances(binA, source)
    n = size(binA, 1);
    d = inf(n, 1);
    d(source) = 0;
    queue = source;
    while ~isempty(queue)
        v = queue(1);
        queue(1) = [];
        nbrs = find(binA(v, :));
        for u = nbrs
            if isinf(d(u))
                d(u) = d(v) + 1;
                queue(end+1) = u; %#ok<AGROW>
            end
        end
    end
end
