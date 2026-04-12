function [adjP, rejected] = bh_fdr(pValues, q)
%BH_FDR Benjamini-Hochberg 1995 step-up FDR adjustment.
%
%   [adjP, rejected] = BH_FDR(pValues, q) returns BH-adjusted p-values
%   (monotonically increasing) and a logical vector indicating which
%   hypotheses pass the false-discovery-rate cut-off Q.
%
% INPUTS:
%   pValues  -  Vector of raw two-sided p-values for an independent
%               family of hypotheses.
%   q        -  Target FDR level (default 0.05).
%
% OUTPUTS:
%   adjP      -  Vector of adjusted p-values, same order as pValues.
%   rejected  -  Logical vector; rejected(i) = adjP(i) <= q.
%
% Reference:
%   Benjamini Y, Hochberg Y. Controlling the false discovery rate: a
%   practical and powerful approach to multiple testing. J R Stat Soc B.
%   1995;57(1):289-300.

    if nargin < 2 || isempty(q); q = 0.05; end
    p = pValues(:);
    m = numel(p);
    if m == 0
        adjP = pValues;
        rejected = false(size(pValues));
        return;
    end

    [pSorted, order] = sort(p, 'ascend');
    ranks = (1:m).';
    adjSorted = pSorted .* (m ./ ranks);

    % Enforce monotonicity (step-up): running min from the top.
    for k = m-1:-1:1
        adjSorted(k) = min(adjSorted(k), adjSorted(k+1));
    end
    adjSorted = min(adjSorted, 1);

    adjP = zeros(size(p));
    adjP(order) = adjSorted;
    adjP = reshape(adjP, size(pValues));
    rejected = adjP <= q;
end
