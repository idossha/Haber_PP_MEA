function result = permutation_test_exact(groupA, groupB, varargin)
%PERMUTATION_TEST_EXACT Exact two-sample permutation test for small n.
%
%   result = PERMUTATION_TEST_EXACT(groupA, groupB) enumerates every
%   possible assignment of the pooled (groupA, groupB) values into
%   groups of size nA and nB, computes a test statistic on each
%   permutation, and returns the exact two-sided p-value. This is the
%   correct test when one or both arms are very small (e.g. n = 3 vs
%   n = 6) and asymptotic tests (t, Mann-Whitney) lose power.
%
%   For n_A = 6 and n_B = 3 there are C(9, 3) = 84 unique assignments,
%   so the test is trivially computable exactly.
%
%   PERMUTATION_TEST_EXACT(groupA, groupB, 'statistic', fcn) uses a
%   custom statistic fcn(a, b) -> scalar. Defaults to the difference of
%   medians median(a) - median(b).
%
% INPUTS:
%   groupA  -  Vector of per-subject values for arm A (e.g. DOI alone).
%   groupB  -  Vector of per-subject values for arm B (e.g. DOI + Ket).
%
% Name-value options:
%   'statistic'  -  Function handle @(a,b) -> scalar. Default:
%                   @(a,b) median(a) - median(b)
%   'tails'      -  'two' (default), 'right', 'left'
%
% OUTPUTS:
%   result  -  Struct with fields:
%       .observed      -  Test statistic on the observed groups
%       .nPermutations -  Number of permutations enumerated
%       .pValue        -  Exact two-sided p-value
%       .pRight        -  P(stat >= observed)
%       .pLeft         -  P(stat <= observed)
%       .permStats     -  Vector of all permutation statistics
%       .tails         -  Tail option used
%
% Reference:
%   Good PI. Permutation, Parametric and Bootstrap Tests of Hypotheses.
%   3rd ed. Springer; 2005. (Also used in Varley et al. 2024, Network
%   Neuroscience, for n = 3 per-condition MEA comparisons.)

    p = inputParser;
    addRequired(p,  'groupA');
    addRequired(p,  'groupB');
    addParameter(p, 'statistic', @(a,b) median(a) - median(b), @(f) isa(f, 'function_handle'));
    addParameter(p, 'tails',     'two', @(s) any(strcmpi(s, {'two','right','left'})));
    parse(p, groupA, groupB, varargin{:});
    opt = p.Results;

    a = groupA(:);
    b = groupB(:);
    a = a(~isnan(a));
    b = b(~isnan(b));
    nA = numel(a);
    nB = numel(b);
    if nA < 1 || nB < 1
        error('permutation_test_exact:Args', 'Each group needs at least one value.');
    end

    pooled = [a; b];
    nTotal = numel(pooled);
    assignmentsA = nchoosek(1:nTotal, nA);     % each row is one permutation
    nPerms       = size(assignmentsA, 1);

    observed = opt.statistic(a, b);

    permStats = zeros(nPerms, 1);
    allIdx    = 1:nTotal;
    for k = 1:nPerms
        idxA = assignmentsA(k, :);
        idxB = setdiff(allIdx, idxA);
        permStats(k) = opt.statistic(pooled(idxA), pooled(idxB));
    end

    pRight = mean(permStats >= observed);
    pLeft  = mean(permStats <= observed);
    switch lower(opt.tails)
        case 'two'
            pValue = 2 * min(pRight, pLeft);
            pValue = min(pValue, 1);
        case 'right'
            pValue = pRight;
        case 'left'
            pValue = pLeft;
    end

    result.observed      = observed;
    result.nPermutations = nPerms;
    result.pValue        = pValue;
    result.pRight        = pRight;
    result.pLeft         = pLeft;
    result.permStats     = permStats;
    result.tails         = lower(opt.tails);
end
