function g = hedges_g_av(baseline, treatment)
%HEDGES_G_AV Small-sample-corrected paired effect size (Lakens 2013).
%
%   g = HEDGES_G_AV(baseline, treatment) computes the paired Hedges'
%   g_av effect size for two within-subject condition vectors of equal
%   length. This is Cohen's d_av (Cumming 2012) with the Hedges small-n
%   bias correction:
%
%       g_av = (mean(t) - mean(b)) / SD_av * (1 - 3 / (4*(n-1) - 1))
%
%   where SD_av = (SD(b) + SD(t)) / 2 (the average of the two condition
%   standard deviations). This is the recommended effect size for paired
%   designs with small n in Lakens 2013.
%
% INPUTS:
%   baseline   -  Vector of baseline values (per-well, per-channel, ...)
%   treatment  -  Matching treatment values (same length).
%
% OUTPUTS:
%   g  -  Scalar Hedges' g_av.
%
% Reference:
%   Lakens D. Calculating and reporting effect sizes to facilitate
%   cumulative science: a practical primer for t-tests and ANOVAs.
%   Front Psychol. 2013;4:863. doi:10.3389/fpsyg.2013.00863

    b = baseline(:);
    t = treatment(:);
    if numel(b) ~= numel(t)
        error('hedges_g_av:Args', 'baseline and treatment must have equal length.');
    end
    ok = ~isnan(b) & ~isnan(t);
    b = b(ok); t = t(ok);
    n = numel(b);
    if n < 2
        g = NaN;
        return;
    end

    sdAv = (std(b) + std(t)) / 2;
    if sdAv == 0
        g = NaN;
        return;
    end

    cohenDav = (mean(t) - mean(b)) / sdAv;
    correction = 1 - 3 / (4 * (n - 1) - 1);
    g = cohenDav * correction;
end
