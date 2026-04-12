function [keep, info] = robust_outlier_filter(yBaseline, yTreatment, varargin)
%ROBUST_OUTLIER_FILTER Drop paired channel observations with extreme values.
%
%   [keep, info] = ROBUST_OUTLIER_FILTER(yBaseline, yTreatment) returns a
%   logical mask KEEP selecting paired channel observations that survive
%   the chosen outlier rule. Default: Tukey 1.5 * IQR **upper** fence on
%   the pooled (baseline; treatment) distribution. Lower fence is never
%   applied because rate metrics are non-negative and small values are
%   legitimate.
%
%   [keep, info] = ROBUST_OUTLIER_FILTER(..., 'mode', 'tukey') — keep only
%       paired observations whose baseline AND treatment are <= the
%       pooled upper fence. Default.
%   [keep, info] = ROBUST_OUTLIER_FILTER(..., 'mode', 'percentile', ...
%                      'upperPercentile', 99) — keep only observations
%       with both values <= the requested percentile of the pooled data.
%   [keep, info] = ROBUST_OUTLIER_FILTER(..., 'mode', 'mean_multiplier', ...
%                      'multiplier', 15) — legacy mode: keep observations
%       whose baseline AND treatment are <= multiplier * mean(condition).
%   [keep, info] = ROBUST_OUTLIER_FILTER(..., 'mode', 'none') — no filter.
%
% INPUTS:
%   yBaseline   -  Baseline values, one per channel.
%   yTreatment  -  Treatment values, one per channel (same length).
%
% Name-value options:
%   'mode'             -  'tukey' (default) | 'percentile' |
%                         'mean_multiplier' | 'none'
%   'upperPercentile'  -  percentile cut for 'percentile' mode (default 99)
%   'multiplier'       -  multiplier for 'mean_multiplier' mode (default 15)
%   'iqrFactor'        -  multiplier on IQR for 'tukey' mode
%                         (default 3.0 = Tukey's "extreme outlier" rule;
%                          1.5 = the standard "mild outlier" rule).
%
% OUTPUTS:
%   keep  -  Logical column vector of length numel(yBaseline). Paired
%            observations to retain.
%   info  -  Struct with fields:
%       .mode              mode used
%       .nTotal            original sample size
%       .nKept             number retained
%       .nDropped          number dropped
%       .upperCutBaseline  upper fence applied to the baseline side
%       .upperCutTreatment upper fence applied to the treatment side
%
% Notes:
%   - The filter is applied symmetrically on the pooled data so that
%     the same cutoff is used across conditions. The cutoff itself is
%     **shared** between baseline and treatment (Mossink 2021 style)
%     rather than computed separately per condition.
%   - For 'mean_multiplier' mode the cutoffs are per-condition because
%     that is how the legacy POC figures defined it.

    p = inputParser;
    addRequired(p,  'yBaseline');
    addRequired(p,  'yTreatment');
    addParameter(p, 'mode',            'tukey', ...
        @(s) any(strcmpi(s, {'tukey','percentile','mean_multiplier','none'})));
    addParameter(p, 'upperPercentile', 99,  @(x) isscalar(x) && x > 50 && x <= 100);
    addParameter(p, 'multiplier',      15,  @(x) isscalar(x) && x > 1);
    addParameter(p, 'iqrFactor',       3.0, @(x) isscalar(x) && x > 0);
    parse(p, yBaseline, yTreatment, varargin{:});
    opt = p.Results;

    yBaseline  = yBaseline(:);
    yTreatment = yTreatment(:);
    n = numel(yBaseline);

    info.mode     = lower(opt.mode);
    info.nTotal   = n;
    info.upperCutBaseline  = NaN;
    info.upperCutTreatment = NaN;

    switch info.mode
        case 'none'
            keep = true(n, 1);

        case 'tukey'
            pooled = [yBaseline; yTreatment];
            pooled = pooled(isfinite(pooled));
            q1 = prctile(pooled, 25);
            q3 = prctile(pooled, 75);
            cut = q3 + opt.iqrFactor * (q3 - q1);
            info.upperCutBaseline  = cut;
            info.upperCutTreatment = cut;
            keep = (yBaseline <= cut) & (yTreatment <= cut);

        case 'percentile'
            pooled = [yBaseline; yTreatment];
            pooled = pooled(isfinite(pooled));
            cut = prctile(pooled, opt.upperPercentile);
            info.upperCutBaseline  = cut;
            info.upperCutTreatment = cut;
            keep = (yBaseline <= cut) & (yTreatment <= cut);

        case 'mean_multiplier'
            cutB = opt.multiplier * mean(yBaseline,  'omitnan');
            cutT = opt.multiplier * mean(yTreatment, 'omitnan');
            info.upperCutBaseline  = cutB;
            info.upperCutTreatment = cutT;
            keep = (yBaseline <= cutB) & (yTreatment <= cutT);
    end

    keep = keep(:);
    info.nKept    = sum(keep);
    info.nDropped = n - info.nKept;
end
