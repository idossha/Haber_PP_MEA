function out = paired_stats(baseline, treatment, varargin)
%PAIRED_STATS All-in-one paired-sample summary for within-well metrics.
%
%   out = PAIRED_STATS(baseline, treatment) computes a canonical set of
%   paired-sample statistics for a single endpoint: median change and
%   percent change, two-level hierarchical-bootstrap 95% CI + p-value,
%   exact Wilcoxon signed-rank p-value (as a non-parametric robustness
%   check), and Hedges' g_av effect size. This is the single function
%   that every figure and table in the paper should call when reporting
%   a paired endpoint, so that reporting is uniform.
%
%   PAIRED_STATS(baseline, treatment, 'nBootstrap', 10000, ...
%        'bootstrapSeed', 1, 'datasetIndex', idx)
%   overrides defaults. If 'datasetIndex' is provided it must be a
%   vector the same length as baseline/treatment assigning each
%   observation to a cluster (typically a recording ID). The
%   hierarchical bootstrap then resamples clusters first, then
%   observations within resampled clusters (Saravanan, Berman & Sober
%   2020). If 'datasetIndex' is empty, a simple one-level bootstrap is
%   used.
%
% INPUTS:
%   baseline   -  Vector of per-observation baseline values.
%   treatment  -  Matching vector of per-observation treatment values.
%
% Name-value options:
%   'nBootstrap'     -  Bootstrap iterations, default 10000.
%   'bootstrapSeed'  -  RNG seed (fixed for reproducibility), default 1.
%   'datasetIndex'   -  (optional) cluster id per observation.
%
% OUTPUTS:
%   out  -  Struct with fields:
%       .n              -  number of valid paired observations
%       .nDatasets      -  number of clusters (or NaN if not clustered)
%       .medianBaseline
%       .medianTreatment
%       .medianDelta            -  treatment - baseline
%       .medianPctChange        -  100 * (t - b) / b (silent-pair excluded)
%       .bootstrap.nIterations
%       .bootstrap.ciDelta       -  95% CI on median delta
%       .bootstrap.ciPctChange   -  95% CI on median percent change
%       .bootstrap.pDelta        -  two-sided percentile p-value
%       .bootstrap.pPctChange    -  two-sided percentile p-value
%       .wilcoxon.W              -  signed-rank statistic
%       .wilcoxon.p              -  exact two-sided p-value
%       .wilcoxon.rRankBiserial  -  rank-biserial effect size
%       .hedgesGav               -  Hedges' g_av (Lakens 2013)
%
% References:
%   Saravanan V, Berman GJ, Sober SJ. Application of the hierarchical
%     bootstrap to multi-level data in neuroscience. NBDT. 2020.
%   Aarts E et al. A solution to dependency. Nat Neurosci. 2014;17:491.
%   Lakens D. Calculating and reporting effect sizes... Front Psychol.
%     2013;4:863.
%
% See also: HEDGES_G_AV, PERMUTATION_TEST_EXACT, BH_FDR.

    p = inputParser;
    addRequired(p,  'baseline');
    addRequired(p,  'treatment');
    addParameter(p, 'nBootstrap',    10000, @(x) isscalar(x) && x > 0);
    addParameter(p, 'bootstrapSeed', 1,     @(x) isscalar(x));
    addParameter(p, 'datasetIndex',  [],    @(x) isempty(x) || isnumeric(x));
    parse(p, baseline, treatment, varargin{:});
    opt = p.Results;

    b = baseline(:);
    t = treatment(:);
    if numel(b) ~= numel(t)
        error('paired_stats:Args', 'baseline and treatment must have equal length.');
    end
    ok = ~isnan(b) & ~isnan(t);
    b = b(ok); t = t(ok);
    if ~isempty(opt.datasetIndex)
        dIdx = opt.datasetIndex(:);
        if numel(dIdx) ~= numel(ok)
            error('paired_stats:Args', ...
                'datasetIndex must match baseline/treatment length.');
        end
        dIdx = dIdx(ok);
    else
        dIdx = [];
    end

    out.n = numel(b);
    if out.n == 0
        out = empty_out();
        return;
    end
    out.nDatasets = iif(isempty(dIdx), NaN, numel(unique(dIdx)));

    out.medianBaseline  = median(b);
    out.medianTreatment = median(t);
    out.medianDelta     = median(t - b);

    % Silent-pair exclusion for percent change (b == 0 & t == 0)
    keepPct = ~(b == 0 & t == 0);
    bP = b(keepPct); tP = t(keepPct);
    pctAll = safe_pct(bP, tP);
    if isempty(pctAll)
        out.medianPctChange = NaN;
    else
        out.medianPctChange = median(pctAll);
    end

    % --- Hierarchical bootstrap ----------------------------------------
    rng(opt.bootstrapSeed, 'twister');
    nIter = opt.nBootstrap;
    bootDelta = nan(nIter, 1);
    bootPct   = nan(nIter, 1);

    if isempty(dIdx)
        for ii = 1:nIter
            idx = randi(out.n, out.n, 1);
            bootDelta(ii) = median(t(idx) - b(idx));
            bBB = b(idx); tBB = t(idx);
            pctBB = safe_pct(bBB, tBB);
            if ~isempty(pctBB)
                bootPct(ii) = median(pctBB);
            end
        end
    else
        clusters = unique(dIdx);
        for ii = 1:nIter
            % Level 1: resample clusters with replacement
            clusSel = clusters(randi(numel(clusters), numel(clusters), 1));
            bB = []; tB = [];
            for cc = 1:numel(clusSel)
                inClust = find(dIdx == clusSel(cc));
                if isempty(inClust); continue; end
                % Level 2: resample observations within cluster
                pick = inClust(randi(numel(inClust), numel(inClust), 1));
                bB = [bB; b(pick)]; %#ok<AGROW>
                tB = [tB; t(pick)]; %#ok<AGROW>
            end
            if isempty(bB); continue; end
            bootDelta(ii) = median(tB - bB);
            pctBB = safe_pct(bB, tB);
            if ~isempty(pctBB)
                bootPct(ii) = median(pctBB);
            end
        end
    end

    out.bootstrap.nIterations = nIter;
    out.bootstrap.ciDelta     = prctile(bootDelta, [2.5 97.5]);
    out.bootstrap.ciPctChange = prctile(bootPct,   [2.5 97.5]);
    out.bootstrap.pDelta      = two_sided_percentile_p(bootDelta);
    out.bootstrap.pPctChange  = two_sided_percentile_p(bootPct);

    % --- Wilcoxon signed-rank (exact for small n, asymptotic otherwise) -
    % signrank + tiedrank live in the Statistics and Machine Learning
    % Toolbox. If that toolbox is not installed we fall back to NaN and
    % add a .wilcoxon.note field; the primary test is the hierarchical
    % bootstrap above, so this is a soft degradation.
    out.wilcoxon.W             = NaN;
    out.wilcoxon.p             = NaN;
    out.wilcoxon.rRankBiserial = NaN;
    out.wilcoxon.note          = '';
    if out.n >= 1 && exist('signrank', 'file') == 2
        try
            try
                [wP, ~, wStats] = signrank(t, b, 'method', 'exact');
            catch
                [wP, ~, wStats] = signrank(t, b);
            end
            out.wilcoxon.W = wStats.signedrank;
            out.wilcoxon.p = wP;
            if exist('tiedrank', 'file') == 2
                out.wilcoxon.rRankBiserial = wilcoxon_r_rank_biserial(t - b);
            else
                out.wilcoxon.note = 'tiedrank unavailable; rank-biserial skipped';
            end
        catch ME
            out.wilcoxon.note = sprintf('signrank failed: %s', ME.message);
        end
    else
        out.wilcoxon.note = 'Statistics Toolbox (signrank) not installed';
    end

    % --- Hedges' g_av ---------------------------------------------------
    out.hedgesGav = hedges_g_av(b, t);
end

% =========================================================================
function out = empty_out()
    out.n = 0;
    out.nDatasets = NaN;
    out.medianBaseline  = NaN;
    out.medianTreatment = NaN;
    out.medianDelta     = NaN;
    out.medianPctChange = NaN;
    out.bootstrap = struct('nIterations', 0, ...
        'ciDelta', [NaN NaN], 'ciPctChange', [NaN NaN], ...
        'pDelta', NaN, 'pPctChange', NaN);
    out.wilcoxon  = struct('W', NaN, 'p', NaN, 'rRankBiserial', NaN);
    out.hedgesGav = NaN;
end

% =========================================================================
function pct = safe_pct(b, t)
    bNZ = b ~= 0;
    pct = nan(size(b));
    pct(bNZ) = 100 * (t(bNZ) - b(bNZ)) ./ b(bNZ);
    % Keep b==0, t>0 as a large positive; drop b==0, t==0 (silent)
    bZtP = (b == 0) & (t > 0);
    pct(bZtP) = 1e6;
    silent = (b == 0) & (t == 0);
    pct(silent) = NaN;
    pct = pct(~isnan(pct));
end

% =========================================================================
function p = two_sided_percentile_p(vec)
    v = vec(~isnan(vec));
    if isempty(v); p = NaN; return; end
    pRight = mean(v >= 0);
    pLeft  = mean(v <= 0);
    p = 2 * min(pRight, pLeft);
    p = min(p, 1);
end

% =========================================================================
function r = wilcoxon_r_rank_biserial(diffs)
% Rank-biserial correlation = (W+ - W-) / (W+ + W-)
    d = diffs(~isnan(diffs) & diffs ~= 0);
    if isempty(d); r = NaN; return; end
    r0 = tiedrank(abs(d));
    wPlus  = sum(r0(d > 0));
    wMinus = sum(r0(d < 0));
    denom = wPlus + wMinus;
    if denom == 0; r = 0; return; end
    r = (wPlus - wMinus) / denom;
end

% =========================================================================
function v = iif(cond, a, b)
    if cond; v = a; else, v = b; end
end
