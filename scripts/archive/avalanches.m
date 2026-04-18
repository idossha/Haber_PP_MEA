function result = avalanches(spikeCells, durationSec, varargin)
%AVALANCHES Neuronal-avalanche statistics (Pasquale 2008, Massobrio 2015).
%
%   result = AVALANCHES(spikeCells, durationSec) detects neuronal
%   avalanches from a pooled-spike timeline binned at the recording's
%   mean inter-event interval (IEI), fits power-law exponents by BOTH
%   least-squares log-log regression (Pasquale 2008 convention) AND
%   Clauset-Shalizi-Newman MLE (Massobrio 2015 convention), classifies
%   the recording as subcritical / critical / supercritical based on
%   distribution shape + GoF, and reports the Sethna scaling residual
%   (Massobrio 2015 Eq. 1-2, p.9) as a criticality-consistency scalar.
%
%   sigma (branching ratio) and DCC are NOT reported: see
%   Tracks/Active/phase1_decisions.md §1 -- neither Pasquale 2008 nor
%   Massobrio 2015 defines these, and the project convention requires
%   every default to be anchored to a page citation in literature/.
%
% INPUTS:
%   spikeCells    -  1 x nCh cell of column-vector spike times (seconds)
%   durationSec   -  scalar recording duration (seconds)
%
% Name-value options (defaults from cfg.avalanches.*):
%   'binMs'             -  override bin width in ms; default: per-recording
%                          mean IEI (Massobrio 2015 p.2, p.8 mean-IEI rule;
%                          Pasquale 2008 p.1357 sweep)
%   'binSweepMs'        -  bin sweep vector (default [0.2 0.4 0.6 0.8 1 2
%                          4 8 16] per Pasquale 2008 p.1357)
%   'runSweep'          -  logical, run the bin sweep for the SI figure
%                          (default false)
%   'tailMassFrac'      -  supercritical tail-mass cutoff fraction of
%                          n_electrodes (default 0.80; Pasquale 2008
%                          Fig. 2/6 upturn rule)
%   'ksPvalueThreshold' -  Clauset KS GoF threshold (default 0.10;
%                          Massobrio 2015 p.5)
%   'alphaSizeTarget'   -  target critical exponent on size distribution
%                          (default -1.5; Pasquale 2008 p.1358)
%   'betaLifetimeTarget'-  target lifetime exponent (default -2.0;
%                          Pasquale 2008 p.1358)
%   'alphaTolerance'    -  |alpha - target| band for "near-critical"
%                          (default 0.15; Pasquale 2008 p.1361 empirical
%                          SD 0.09-0.13)
%
% OUTPUTS:
%   result.binMs                  -  bin used for the primary fit (ms)
%   result.meanIeiMs              -  computed mean IEI across all spikes (ms)
%   result.nAvalanches            -  total count at the primary bin
%   result.sizeCounts             -  histogram of sizes (integer bins 1..nCh)
%   result.lifetimeCounts         -  histogram of lifetimes (integer bins)
%   result.alphaLS                -  LS log-log size exponent (Pasquale p.1357)
%   result.betaLS                 -  LS log-log lifetime exponent
%   result.alphaLS_rmse           -  RMSE of LS size fit in log-log space
%   result.betaLS_rmse            -  RMSE of LS lifetime fit in log-log space
%   result.alphaMLE               -  Clauset MLE size exponent (positive-valued)
%   result.alphaMLE_xmin          -  MLE x_min from KS minimisation
%   result.alphaMLE_ksp           -  KS p-value from 100-iter semi-parametric bootstrap
%   result.lrtPreferredModel      -  'power_law' | 'exponential' |
%                                    'truncated_power_law' | 'lognormal'
%   result.tailMassFrac           -  fraction of avalanches with size >=
%                                    tailMassFrac*nCh
%   result.sethnaResidual         -  |(tau-1)/(alpha-1) - 1/(sigma*nu*z)_fit|
%                                    (Massobrio 2015 Eq. 2 p.9)
%   result.classification         -  'subcritical' | 'critical' |
%                                    'supercritical' | 'indeterminate'
%   result.binSweep               -  struct array (only populated if
%                                    'runSweep', true); fields: binMs,
%                                    alphaLS, betaLS
%
% References:
%   Pasquale V, Massobrio P, Bologna LL, Chiappalone M, Martinoia S.
%     Self-organization and neuronal avalanches in networks of dissociated
%     cortical neurons. Neuroscience 2008;153:1354. Methods p.1357
%     (avalanche defn, LS fit + drop-x1 + drop-<1%max). Results p.1358,
%     p.1361 (critical exponents, empirical spread). Discussion p.1367
%     (sub/critical/super shape rules).
%   Massobrio P, Pasquale V, Martinoia S. Self-organized criticality in
%     cortical assemblies. Sci Rep 2015;5:10578. p.2 (active-bin rule,
%     sub/supercritical shape rules). p.5 (KS p-value >= 0.1 threshold).
%     p.9 Eq. 1-2 (Sethna scaling relation). p.14 (Clauset MLE + LRT
%     alternatives set).
%   Clauset A, Shalizi CR, Newman MEJ. Power-law distributions in
%     empirical data. SIAM Rev 2009;51:661. (Discrete MLE and KS
%     semi-parametric bootstrap recipe cited inline in Massobrio 2015 p.14.)
%
% See also: TRANSFER_ENTROPY_D1, BURST_ONSET_COINCIDENCE, PAIRED_STATS.

    %%% ----- parse inputs --------------------------------------------------
    cfg = project_config();
    defBinSweep = pick_default(cfg, 'avalanches', 'bin_sweep_ms', ...
        [0.2 0.4 0.6 0.8 1 2 4 8 16]);
    defTailMass = pick_default(cfg, 'avalanches', 'tail_mass_frac', 0.80);
    defKSthr    = pick_default(cfg, 'avalanches', 'ks_pvalue_threshold', 0.10);
    defAlphaTgt = pick_default(cfg, 'avalanches', 'alpha_size_target', -1.5);
    defBetaTgt  = pick_default(cfg, 'avalanches', 'beta_lifetime_target', -2.0);
    defAlphaTol = pick_default(cfg, 'avalanches', 'alpha_tolerance', 0.15);

    p = inputParser;
    addRequired(p,  'spikeCells',   @iscell);
    addRequired(p,  'durationSec',  @(x) isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'binMs',              [],          @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'binSweepMs',         defBinSweep, @(x) isnumeric(x) && isvector(x) && all(x > 0));
    addParameter(p, 'runSweep',           false,       @(x) isscalar(x) && islogical(x));
    addParameter(p, 'tailMassFrac',       defTailMass, @(x) isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'ksPvalueThreshold',  defKSthr,    @(x) isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'alphaSizeTarget',    defAlphaTgt, @(x) isscalar(x));
    addParameter(p, 'betaLifetimeTarget', defBetaTgt,  @(x) isscalar(x));
    addParameter(p, 'alphaTolerance',     defAlphaTol, @(x) isscalar(x) && x > 0);
    parse(p, spikeCells, durationSec, varargin{:});
    opt = p.Results;

    nCh = numel(spikeCells);
    result = empty_avalanche_result(nCh);

    %%% ----- pool spikes and compute mean IEI -----------------------------
    % Pasquale 2008 p.1357 / Massobrio 2015 p.2: active bin = at least one
    % spike anywhere on the array. Build a pooled, sorted event list first
    % so the mean IEI is computed over ALL spikes (Massobrio 2015 p.2,
    % p.8 mean-IEI rule).
    pooled = [];
    for c = 1:nCh
        st = spikeCells{c};
        if ~isempty(st); pooled = [pooled; st(:)]; end %#ok<AGROW>
    end
    pooled = sort(pooled);
    if numel(pooled) < 20
        return;
    end

    allIEIs = diff(pooled);
    meanIeiMs = mean(allIEIs) * 1000;
    result.meanIeiMs = meanIeiMs;

    %%% ----- bin selection -------------------------------------------------
    % If the caller supplied a bin, use it verbatim. Otherwise snap the
    % mean IEI to the nearest bin in binSweepMs and floor at 0.2 ms
    % (Pasquale 2008 p.1360 reports IEI ~0.3 ms for dense cultures, and
    % 0.2 ms is the smallest bin in the canonical sweep).
    if isempty(opt.binMs)
        targetMs = max(0.2, meanIeiMs);
        [~, bi] = min(abs(opt.binSweepMs - targetMs));
        primaryBinMs = opt.binSweepMs(bi);
    else
        primaryBinMs = opt.binMs;
    end
    result.binMs = primaryBinMs;

    %%% ----- primary avalanche detection + statistics --------------------
    primary = avalanche_stats(spikeCells, durationSec, pooled, primaryBinMs, nCh, opt);

    % Copy primary fields to the top-level result
    result.nAvalanches       = primary.nAvalanches;
    result.sizeCounts        = primary.sizeCounts;
    result.lifetimeCounts    = primary.lifetimeCounts;
    result.alphaLS           = primary.alphaLS;
    result.betaLS            = primary.betaLS;
    result.alphaLS_rmse      = primary.alphaLS_rmse;
    result.betaLS_rmse       = primary.betaLS_rmse;
    result.alphaMLE          = primary.alphaMLE;
    result.alphaMLE_xmin     = primary.alphaMLE_xmin;
    result.alphaMLE_ksp      = primary.alphaMLE_ksp;
    result.lrtPreferredModel = primary.lrtPreferredModel;
    result.tailMassFrac      = primary.tailMassFrac;
    result.sethnaResidual    = primary.sethnaResidual;

    %%% ----- classification tree ------------------------------------------
    result.classification = classify_regime(primary, opt);

    %%% ----- optional bin sweep for SI §S2 --------------------------------
    if opt.runSweep
        sweep = struct('binMs', {}, 'alphaLS', {}, 'betaLS', {});
        for k = 1:numel(opt.binSweepMs)
            bms = opt.binSweepMs(k);
            try
                stmp = avalanche_stats(spikeCells, durationSec, pooled, bms, nCh, opt);
                sweep(end + 1) = struct( ...
                    'binMs', bms, 'alphaLS', stmp.alphaLS, ...
                    'betaLS', stmp.betaLS); %#ok<AGROW>
            catch
                sweep(end + 1) = struct( ...
                    'binMs', bms, 'alphaLS', NaN, 'betaLS', NaN); %#ok<AGROW>
            end
        end
        result.binSweep = sweep;
    end
end

% =========================================================================
function out = avalanche_stats(spikeCells, durationSec, pooled, binMs, nCh, opt)
%AVALANCHE_STATS Run a single-bin avalanche detection + fitting pipeline.

    out = empty_stats_block();
    binSec = binMs / 1000;
    edges = 0:binSec:durationSec;
    nBins = numel(edges) - 1;
    if nBins < 10; return; end

    %%% ----- detect active bins -------------------------------------------
    % Pasquale 2008 p.1357 / Massobrio 2015 p.2: bin is active iff at
    % least one spike anywhere on the array falls in it.
    activePooled = histcounts(pooled, edges) > 0;

    % Per-channel binary matrix for later size computation.
    chanBins = false(nCh, nBins);
    for c = 1:nCh
        st = spikeCells{c};
        if isempty(st); continue; end
        chanBins(c, :) = histcounts(st, edges) > 0;
    end

    %%% ----- run-length-encode active bins into avalanches ---------------
    % Each contiguous run of true bins is one avalanche. Drop runs that
    % touch the recording boundary per Pasquale 2008 p.1357 flanking-
    % silent-bin rule.
    d = diff([false, activePooled, false]);
    startBins = find(d == 1);
    endBins   = find(d == -1) - 1;

    % Boundary-truncation: drop runs starting at bin 1 or ending at nBins.
    keep = true(size(startBins));
    if ~isempty(startBins)
        if startBins(1) == 1;        keep(1)   = false; end
        if endBins(end) == nBins;    keep(end) = false; end
    end
    startBins = startBins(keep);
    endBins   = endBins(keep);

    if isempty(startBins); return; end

    %%% ----- size and lifetime per avalanche ------------------------------
    % Pasquale 2008 definition 2 (p.1357, p.1361): size = #unique electrodes
    % active at least once during the avalanche.
    nAv = numel(startBins);
    sizes = zeros(nAv, 1);
    lifetimes = zeros(nAv, 1);
    for k = 1:nAv
        s = startBins(k);
        e = endBins(k);
        % Any channel firing at least once in bins [s..e]
        sub = any(chanBins(:, s:e), 2);
        sizes(k) = sum(sub);
        lifetimes(k) = e - s + 1;
    end

    % Drop avalanches with zero size (shouldn't happen, but be defensive).
    ok = sizes >= 1;
    sizes = sizes(ok);
    lifetimes = lifetimes(ok);
    if isempty(sizes); return; end

    out.nAvalanches = numel(sizes);

    %%% ----- histograms ---------------------------------------------------
    sizeEdges = 0.5:1:(nCh + 0.5);
    out.sizeCounts = histcounts(sizes, sizeEdges);
    maxLife = max(lifetimes);
    lifeEdges = 0.5:1:(maxLife + 0.5);
    out.lifetimeCounts = histcounts(lifetimes, lifeEdges);

    %%% ----- LS log-log fit (Pasquale 2008 p.1357) ------------------------
    [out.alphaLS, out.alphaLS_rmse] = fit_powerlaw_ls(out.sizeCounts);
    [out.betaLS,  out.betaLS_rmse]  = fit_powerlaw_ls(out.lifetimeCounts);

    %%% ----- Clauset MLE (Massobrio 2015 p.14) ----------------------------
    try
        [mleAlpha, mleXmin, mleKsp] = fit_powerlaw_mle(sizes);
        out.alphaMLE      = mleAlpha;
        out.alphaMLE_xmin = mleXmin;
        out.alphaMLE_ksp  = mleKsp;
    catch
        out.alphaMLE      = NaN;
        out.alphaMLE_xmin = NaN;
        out.alphaMLE_ksp  = NaN;
    end

    %%% ----- LRT vs alternatives (Massobrio 2015 p.14) --------------------
    try
        out.lrtPreferredModel = lrt_preferred_model(sizes, out.alphaMLE, out.alphaMLE_xmin);
    catch
        out.lrtPreferredModel = 'power_law';
    end

    %%% ----- tail mass at system size (Pasquale 2008 Fig. 2/6) -----------
    out.tailMassFrac = sum(sizes >= opt.tailMassFrac * nCh) / numel(sizes);

    %%% ----- Sethna scaling residual (Massobrio 2015 Eq. 1-2, p.9) -------
    out.sethnaResidual = sethna_residual(sizes, lifetimes, out.alphaLS, out.betaLS);
end

% =========================================================================
function [alpha, rmse] = fit_powerlaw_ls(counts)
%FIT_POWERLAW_LS Least-squares log-log regression (Pasquale 2008 p.1357).
%
%   Excludes bin 1 (unit-size / unit-duration) and bins whose count is
%   below 1% of max(count). Returns NaN-safe fit.

    alpha = NaN; rmse = NaN;
    if isempty(counts) || sum(counts) == 0; return; end

    x = 1:numel(counts);
    c = counts(:).';

    % Drop bin 1 (Pasquale 2008 p.1357).
    keep = x >= 2;
    % Drop bins with count < 1% of max (Pasquale 2008 p.1357).
    keep = keep & (c >= 0.01 * max(c));
    % Drop zero-count bins so log is defined.
    keep = keep & (c > 0);

    if sum(keep) < 3; return; end

    lx = log10(x(keep));
    ly = log10(c(keep));
    P  = polyfit(lx, ly, 1);
    alpha = P(1);
    residuals = ly - polyval(P, lx);
    rmse = sqrt(mean(residuals .^ 2));
end

% =========================================================================
function [alpha, xminBest, ksp] = fit_powerlaw_mle(sizes)
%FIT_POWERLAW_MLE Clauset-Shalizi-Newman discrete power-law MLE.
%
%   Returns positive alpha (Clauset convention alpha > 1). KS p-value is
%   estimated via 100-iteration semi-parametric bootstrap.

    alpha = NaN; xminBest = NaN; ksp = NaN;
    s = sizes(:);
    s = s(s >= 1);
    if numel(s) < 20; return; end

    candidates = unique(s);
    candidates = candidates(candidates <= max(s) - 2);   % need tail to fit
    if isempty(candidates); candidates = unique(s); end

    bestD = Inf;
    bestAlpha = NaN;
    bestXmin = NaN;
    for k = 1:numel(candidates)
        xm = candidates(k);
        tail = s(s >= xm);
        n = numel(tail);
        if n < 10; continue; end
        % Continuous-approx MLE (Clauset 2009 eq 3.7) — fine for
        % integer-valued avalanche sizes above a few counts.
        denom = sum(log(tail ./ (xm - 0.5)));
        if denom <= 0; continue; end
        aHat = 1 + n / denom;
        if ~isfinite(aHat) || aHat <= 1; continue; end

        % KS distance on the tail CDF
        D = ks_distance_discrete(tail, aHat, xm);
        if D < bestD
            bestD = D;
            bestAlpha = aHat;
            bestXmin  = xm;
        end
    end
    if ~isfinite(bestAlpha); return; end

    alpha = bestAlpha;
    xminBest = bestXmin;

    %%% ----- semi-parametric bootstrap KS p-value (100 iter) --------------
    nBoot = 100;
    dBoot = nan(nBoot, 1);
    nData = numel(s);
    tailFrac = sum(s >= bestXmin) / nData;
    rng(1, 'twister');
    for ii = 1:nBoot
        synth = zeros(nData, 1);
        isTail = rand(nData, 1) < tailFrac;
        nTail = sum(isTail);
        nBulk = nData - nTail;
        % Bulk: resample from observed bulk
        bulk = s(s < bestXmin);
        if ~isempty(bulk) && nBulk > 0
            synth(~isTail) = bulk(randi(numel(bulk), nBulk, 1));
        end
        % Tail: sample from fitted power law via inverse transform
        if nTail > 0
            synth(isTail) = sample_discrete_powerlaw(nTail, bestAlpha, bestXmin);
        end
        % Refit x_min-aware on synthetic tail (use same xmin for speed)
        tailB = synth(synth >= bestXmin);
        if numel(tailB) < 5
            dBoot(ii) = NaN;
            continue;
        end
        denom = sum(log(tailB ./ (bestXmin - 0.5)));
        if denom <= 0; dBoot(ii) = NaN; continue; end
        aB = 1 + numel(tailB) / denom;
        if ~isfinite(aB) || aB <= 1; dBoot(ii) = NaN; continue; end
        dBoot(ii) = ks_distance_discrete(tailB, aB, bestXmin);
    end
    dBoot = dBoot(~isnan(dBoot));
    if isempty(dBoot)
        ksp = NaN;
    else
        ksp = mean(dBoot >= bestD);
    end
end

% =========================================================================
function D = ks_distance_discrete(tail, alpha, xmin)
%KS_DISTANCE_DISCRETE Max deviation between empirical and fitted CDFs.
    u = sort(unique(tail(:)));
    nTail = numel(tail);
    % Empirical CDF at unique points
    empCDF = arrayfun(@(x) sum(tail <= x) / nTail, u);
    % Fitted CDF: continuous approximation F(x) = 1 - (x/xmin)^(1-alpha)
    fitCDF = 1 - (u ./ xmin) .^ (1 - alpha);
    D = max(abs(empCDF - fitCDF));
end

% =========================================================================
function samples = sample_discrete_powerlaw(n, alpha, xmin)
%SAMPLE_DISCRETE_POWERLAW Inverse-transform sampling via continuous approx.
    u = rand(n, 1);
    samples = round((xmin - 0.5) .* (1 - u) .^ (1 / (1 - alpha)) + 0.5);
    samples(samples < xmin) = xmin;
end

% =========================================================================
function preferred = lrt_preferred_model(sizes, alphaMLE, xmin)
%LRT_PREFERRED_MODEL Vuong-style LRT vs alternatives (Massobrio 2015 p.14).
%
%   Alternatives: exponential, truncated power-law, log-normal. Returns
%   the model with the highest total log-likelihood on the tail.

    preferred = 'power_law';
    s = sizes(:);
    s = s(s >= 1);
    if isempty(s) || ~isfinite(alphaMLE) || ~isfinite(xmin)
        return;
    end
    tail = s(s >= xmin);
    if numel(tail) < 10; return; end

    % Log-likelihoods (per-point, summed)
    llPL  = ll_powerlaw(tail, alphaMLE, xmin);
    llExp = ll_exponential(tail, xmin);
    llTPL = ll_trunc_powerlaw(tail, xmin);
    llLN  = ll_lognormal(tail, xmin);

    lls = [llPL, llExp, llTPL, llLN];
    names = {'power_law', 'exponential', 'truncated_power_law', 'lognormal'};
    [~, idx] = max(lls);
    preferred = names{idx};
end

% =========================================================================
function ll = ll_powerlaw(tail, alpha, xmin)
    % Normalised PL PDF: p(x) = (alpha-1) * xmin^(alpha-1) * x^(-alpha)
    % log p(x) = log(alpha-1) + (alpha-1)*log(xmin) - alpha*log(x)
    if alpha <= 1
        ll = -Inf; return;
    end
    ll = sum(log(alpha - 1) + (alpha - 1) * log(xmin) - alpha * log(tail));
    if ~isfinite(ll); ll = -Inf; end
end

% =========================================================================
function ll = ll_exponential(tail, xmin)
    % MLE for shifted exponential: lambda = 1 / (mean(tail) - xmin + 1)
    m = mean(tail);
    lambda = 1 / max(m - xmin + 1, 1e-6);
    ll = sum(log(lambda) - lambda * (tail - xmin));
    if ~isfinite(ll); ll = -Inf; end
end

% =========================================================================
function ll = ll_trunc_powerlaw(tail, xmin)
    % Truncated PL: exp cut-off. Fit with a coarse grid search over
    % (alpha, lambda) — sufficient for an LRT comparison.
    alphaGrid  = linspace(1.1, 3.5, 10);
    lambdaGrid = logspace(-4, -1, 8);
    bestLL = -Inf;
    for ai = 1:numel(alphaGrid)
        for li = 1:numel(lambdaGrid)
            a = alphaGrid(ai);
            lam = lambdaGrid(li);
            % Unnormalised logpdf: -a*log(x) - lam*x
            logp = -a * log(tail) - lam * tail;
            % Approximate normaliser via sum over x in [xmin, max(tail)]
            xs = xmin:max(tail);
            Z = sum(xs .^ (-a) .* exp(-lam * xs));
            if Z <= 0 || ~isfinite(Z); continue; end
            ll = sum(logp) - numel(tail) * log(Z);
            if ll > bestLL; bestLL = ll; end
        end
    end
    ll = bestLL;
    if ~isfinite(ll); ll = -Inf; end
end

% =========================================================================
function ll = ll_lognormal(tail, xmin) %#ok<INUSD>
    lx = log(tail);
    mu = mean(lx);
    sigma = std(lx);
    if sigma <= 0 || ~isfinite(sigma); ll = -Inf; return; end
    ll = sum(-log(tail .* sigma * sqrt(2*pi)) - (lx - mu).^2 / (2 * sigma^2));
    if ~isfinite(ll); ll = -Inf; end
end

% =========================================================================
function residual = sethna_residual(sizes, lifetimes, alphaLS, betaLS)
%SETHNA_RESIDUAL Massobrio 2015 Eq. 2 (p.9) criticality-consistency check.
%
%   residual = | (tau - 1)/(alpha - 1) - 1/(sigma*nu*z)_fit |
%
%   with alpha = |alphaLS|, tau = |betaLS|, and 1/(sigma*nu*z) from the
%   slope of log<s>(T) vs log T (Massobrio Fig. 5c,f).

    residual = NaN;
    if ~isfinite(alphaLS) || ~isfinite(betaLS); return; end
    a = abs(alphaLS);
    tau = abs(betaLS);
    if a <= 1 || tau <= 1; return; end

    % Bin <s>(T) by lifetime
    Tunique = unique(lifetimes);
    meanS = nan(numel(Tunique), 1);
    for k = 1:numel(Tunique)
        meanS(k) = mean(sizes(lifetimes == Tunique(k)));
    end
    ok = meanS > 0 & Tunique > 0;
    if sum(ok) < 3; return; end

    lT = log(Tunique(ok));
    lS = log(meanS(ok));
    P = polyfit(lT, lS, 1);
    invSNZ_fit = P(1);
    if ~isfinite(invSNZ_fit); return; end

    invSNZ_pred = (tau - 1) / (a - 1);
    residual = abs(invSNZ_pred - invSNZ_fit);
end

% =========================================================================
function label = classify_regime(primary, opt)
%CLASSIFY_REGIME Distribution-shape classification tree.
%
%   Tree from avalanche_method_spec.md §10 / Tracks/Active/phase1_decisions.md §1:
%
%   - critical      = KS p >= 0.10 AND |alphaLS - 1.5| < 0.15 AND LRT == power_law
%   - subcritical   = LRT == exponential OR |alphaLS| > 1.9
%   - supercritical = tailMassFrac > 0.05 OR LRT == truncated_power_law
%   - else          = indeterminate

    label = 'indeterminate';

    alphaAbs = abs(primary.alphaLS);

    criticalHit = isfinite(primary.alphaMLE_ksp) ...
        && primary.alphaMLE_ksp >= opt.ksPvalueThreshold ...
        && isfinite(primary.alphaLS) ...
        && abs(alphaAbs - abs(opt.alphaSizeTarget)) < opt.alphaTolerance ...
        && strcmp(primary.lrtPreferredModel, 'power_law');
    if criticalHit
        label = 'critical';
        return;
    end

    subHit = strcmp(primary.lrtPreferredModel, 'exponential') ...
        || (isfinite(primary.alphaLS) && alphaAbs > 1.9);
    if subHit
        label = 'subcritical';
        return;
    end

    superHit = (isfinite(primary.tailMassFrac) && primary.tailMassFrac > 0.05) ...
        || strcmp(primary.lrtPreferredModel, 'truncated_power_law');
    if superHit
        label = 'supercritical';
        return;
    end
end

% =========================================================================
function r = empty_avalanche_result(nCh)
    r.binMs               = NaN;
    r.meanIeiMs           = NaN;
    r.nAvalanches         = 0;
    r.sizeCounts          = zeros(1, nCh);
    r.lifetimeCounts      = zeros(1, 0);
    r.alphaLS             = NaN;
    r.betaLS              = NaN;
    r.alphaLS_rmse        = NaN;
    r.betaLS_rmse         = NaN;
    r.alphaMLE            = NaN;
    r.alphaMLE_xmin       = NaN;
    r.alphaMLE_ksp        = NaN;
    r.lrtPreferredModel   = 'power_law';
    r.tailMassFrac        = 0;
    r.sethnaResidual      = NaN;
    r.classification      = 'indeterminate';
    r.binSweep            = struct('binMs', {}, 'alphaLS', {}, 'betaLS', {});
end

% =========================================================================
function r = empty_stats_block()
    r.nAvalanches         = 0;
    r.sizeCounts          = zeros(1, 0);
    r.lifetimeCounts      = zeros(1, 0);
    r.alphaLS             = NaN;
    r.betaLS              = NaN;
    r.alphaLS_rmse        = NaN;
    r.betaLS_rmse         = NaN;
    r.alphaMLE            = NaN;
    r.alphaMLE_xmin       = NaN;
    r.alphaMLE_ksp        = NaN;
    r.lrtPreferredModel   = 'power_law';
    r.tailMassFrac        = 0;
    r.sethnaResidual      = NaN;
end

% =========================================================================
function v = pick_default(cfg, section, field, fallback)
    if isfield(cfg, section) && isfield(cfg.(section), field)
        v = cfg.(section).(field);
    else
        v = fallback;
    end
end
