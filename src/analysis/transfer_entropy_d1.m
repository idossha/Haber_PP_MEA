function result = transfer_entropy_d1(spikeCells, durationSec, varargin)
%TRANSFER_ENTROPY_D1 Delay-1 pairwise transfer entropy on binned spike trains.
%
%   result = TRANSFER_ENTROPY_D1(spikeCells, durationSec) computes the
%   pairwise delay-1 transfer entropy matrix TE(i->j) on binary-binned
%   spike trains from a single recording, then returns a summary struct
%   containing the per-recording scalars that feed into paired_stats.
%
%   TE(X->Y) = sum p(y_{t+1}, y_t, x_t) log2 [ p(y_{t+1}|y_t,x_t) / p(y_{t+1}|y_t) ]
%
%   The function first discretises each channel's spike train to binary
%   occupancy at a 1 ms bin (Ito 2011 p.4), gates channels by minimum
%   firing rate (default 5 Hz), then computes D1TE on ALL eligible
%   channels (no top-K subsampling, unlike the tmp/ prototype). The
%   meanTE scalar is the mean over the top decile of off-diagonal entries
%   and the asymmetry scalar is the mean |TE_{ij} - TE_{ji}| normalised
%   by <TE> over those same top-decile ordered pairs.
%
% INPUTS:
%   spikeCells   -  1 x nCh cell of column-vector spike times (seconds)
%   durationSec  -  scalar recording duration (seconds)
%
% Name-value options (defaults from project_config().te):
%   'binMs'         -  bin width in ms (default cfg.te.bin_ms = 1 ms;
%                       Ito 2011 p.4 canonical D1TE bin)
%   'minRateHz'     -  active-channel eligibility threshold in Hz
%                       (default cfg.te.min_rate_hz = 5; Ito 2011
%                       active-channel convention)
%   'topEdgesFrac'  -  top-decile fraction for meanTE / asymmetry
%                       summary (default cfg.te.top_edges_frac = 0.10)
%
% OUTPUTS:
%   result.TE             -  nCh x nCh TE matrix in bits (non-active
%                            channels retained as zero rows/cols so the
%                            caller can index by the original channel id)
%   result.meanTE         -  scalar, mean of top-decile off-diagonal
%                            entries (bits)
%   result.asymmetry      -  scalar, <|TE_{ij} - TE_{ji}|> / <TE> taken
%                            over the top-decile ordered pairs
%   result.activeChannels -  indices of channels passing the minRateHz
%                            eligibility gate
%   result.nActive        -  numel(activeChannels)
%   result.binMs          -  echoed bin width
%   result.minRateHz      -  echoed rate gate
%   result.topEdgesFrac   -  echoed top-decile fraction
%
% References:
%   Ito S, Hansen ME, Heiland R, Lumsdaine A, Litke AM, Beggs JM.
%     Extending transfer entropy improves identification of effective
%     connectivity in a spiking cortical network model. PLoS ONE
%     2011;6:e27431. (p.4 D1TE definition, p.5 1 ms bin convention)
%   Varley TF, Pope M, Grazzi LF, Pittman-Polletta B, Sporns O.
%     Information processing dynamics in neural networks of macroscopic
%     brain organoids. J Neural Eng 2024. (§2.3 mTE analog).
%
% See also: CONNECTIVITY_XCORR, BURST_ONSET_COINCIDENCE, SPIKE_ENTROPY_SHANNON.

    %%% ----- parse inputs --------------------------------------------------
    cfg = project_config();
    defBinMs         = pick_default(cfg, 'te', 'bin_ms',         1);
    defMinRateHz     = pick_default(cfg, 'te', 'min_rate_hz',    5);
    defTopEdgesFrac  = pick_default(cfg, 'te', 'top_edges_frac', 0.10);

    p = inputParser;
    addRequired(p,  'spikeCells',  @iscell);
    addRequired(p,  'durationSec', @(x) isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'binMs',        defBinMs,        @(x) isscalar(x) && x > 0);
    addParameter(p, 'minRateHz',    defMinRateHz,    @(x) isscalar(x) && x >= 0);
    addParameter(p, 'topEdgesFrac', defTopEdgesFrac, @(x) isscalar(x) && x > 0 && x <= 1);
    parse(p, spikeCells, durationSec, varargin{:});
    opt = p.Results;

    nCh = numel(spikeCells);
    result = empty_te_result(nCh, opt);

    %%% ----- bin to binary occupancy --------------------------------------
    % Ito 2011 p.4: 1 ms bins, binary occupancy (>=1 spike -> 1).
    edges = 0:(opt.binMs/1000):durationSec;
    nBins = numel(edges) - 1;
    if nBins < 10
        return;
    end

    M = false(nCh, nBins);
    for c = 1:nCh
        st = spikeCells{c};
        if isempty(st); continue; end
        M(c, :) = histcounts(st, edges) > 0;
    end

    %%% ----- active-channel gate ------------------------------------------
    % Ito 2011 uses a minimum firing-rate gate to exclude essentially-silent
    % channels (they contribute no information and dominate the top-decile
    % mean with numerical zeros).
    counts = cellfun(@numel, spikeCells(:));
    rates  = counts / durationSec;
    active = find(rates >= opt.minRateHz);
    if numel(active) < 4
        return;
    end

    %%% ----- D1TE on all active pairs -------------------------------------
    TE = zeros(nCh, nCh);
    for a = 1:numel(active)
        ia = active(a);
        xa = M(ia, :);
        for b = 1:numel(active)
            if a == b; continue; end
            ib = active(b);
            TE(ia, ib) = d1te(xa, M(ib, :));
        end
    end

    %%% ----- top-decile meanTE and asymmetry ------------------------------
    % Off-diagonal entries restricted to active x active sub-block.
    subTE = TE(active, active);
    nActive = numel(active);
    nPairs  = nActive * (nActive - 1);
    if nPairs < 10
        result.TE = TE;
        result.activeChannels = active(:).';
        result.nActive = nActive;
        return;
    end

    % Threshold picks the top `topEdgesFrac` entries of the off-diagonal
    % values. Same convention as tmp/NEW_ANGLES_REPORT.md §1 summary.
    mask = ~eye(nActive, 'logical');
    vals = subTE(mask);
    vals = sort(vals, 'descend');
    topN = max(5, round(opt.topEdgesFrac * nPairs));
    topN = min(topN, numel(vals));
    thr  = vals(topN);
    topMask = subTE >= thr & mask;

    % meanTE = mean over top-decile entries (positive only, to avoid
    % numerical negatives dragging the mean).
    topVals = subTE(topMask);
    topVals = topVals(topVals > 0);
    if isempty(topVals)
        meanTE = 0;
    else
        meanTE = mean(topVals);
    end

    % asymmetry = mean over unique ordered pairs (i,j) in the top mask of
    % |TE(i,j) - TE(j,i)| divided by the mean of TE over those same pairs.
    [I, J] = find(topMask);
    asymVals = nan(numel(I), 1);
    seen = false(nActive, nActive);
    k = 0;
    for kk = 1:numel(I)
        ii = I(kk); jj = J(kk);
        if seen(ii, jj) || seen(jj, ii); continue; end
        k = k + 1;
        asymVals(k) = abs(subTE(ii, jj) - subTE(jj, ii));
        seen(ii, jj) = true; seen(jj, ii) = true;
    end
    asymVals = asymVals(1:k);
    if isempty(asymVals) || meanTE <= 0
        asymmetry = 0;
    else
        asymmetry = mean(asymVals) / meanTE;
    end

    %%% ----- pack ----------------------------------------------------------
    result.TE              = TE;
    result.meanTE          = meanTE;
    result.asymmetry       = asymmetry;
    result.activeChannels  = active(:).';
    result.nActive         = nActive;
    result.binMs           = opt.binMs;
    result.minRateHz       = opt.minRateHz;
    result.topEdgesFrac    = opt.topEdgesFrac;
end

% =========================================================================
function r = empty_te_result(nCh, opt)
    r.TE              = zeros(nCh, nCh);
    r.meanTE          = 0;
    r.asymmetry       = 0;
    r.activeChannels  = zeros(1, 0);
    r.nActive         = 0;
    r.binMs           = opt.binMs;
    r.minRateHz       = opt.minRateHz;
    r.topEdgesFrac    = opt.topEdgesFrac;
end

% =========================================================================
function te = d1te(x, y)
%D1TE Delay-1 transfer entropy X->Y on binary sequences (bits).
%
% Ported verbatim from tmp/m/tmp_angle_te.m lines 133-157. Uses base-2
% log so units are bits. Non-finite or numerically-negative values are
% clamped to zero.
    y1 = y(2:end); y0 = y(1:end-1); x0 = x(1:end-1);
    idx = 1 + y1 + 2*y0 + 4*x0;   % 0..7 triplet index
    cnt = accumarray(idx(:), 1, [8 1]);
    N = sum(cnt);
    if N == 0; te = 0; return; end
    p = cnt / N;
    p_y0x0 = accumarray(1 + y0(:) + 2*x0(:), 1, [4 1]) / N;
    p_y1y0 = accumarray(1 + y1(:) + 2*y0(:), 1, [4 1]) / N;
    p_y0   = accumarray(1 + y0(:),            1, [2 1]) / N;
    te = 0;
    for a = 1:8
        if p(a) == 0; continue; end
        [y1v, y0v, x0v] = decode3(a - 1);
        p_joint = p(a);
        p_condA = p_joint / max(p_y0x0(1 + y0v + 2*x0v), eps);
        p_condB = p_y1y0(1 + y1v + 2*y0v) / max(p_y0(1 + y0v), eps);
        if p_condA > 0 && p_condB > 0
            te = te + p_joint * log2(p_condA / p_condB);
        end
    end
    if ~isfinite(te) || te < 0; te = max(te, 0); end
end

% =========================================================================
function [y1, y0, x0] = decode3(v)
    y1 = bitand(v, 1);
    y0 = bitand(bitshift(v, -1), 1);
    x0 = bitand(bitshift(v, -2), 1);
end

% =========================================================================
function v = pick_default(cfg, section, field, fallback)
    if isfield(cfg, section) && isfield(cfg.(section), field)
        v = cfg.(section).(field);
    else
        v = fallback;
    end
end
