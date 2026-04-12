function result = burst_onset_coincidence(burstStartTimes, varargin)
%BURST_ONSET_COINCIDENCE Pairwise cross-channel burst-onset simultaneity.
%
%   result = BURST_ONSET_COINCIDENCE(burstStartTimes) computes the mean
%   pairwise burst-onset coincidence across channels. For each channel
%   pair (i, j) with burst start-time sets B_i, B_j:
%
%       C_ij = #{ (b_i, b_j) : |b_i - b_j| <= windowMs } / sqrt(|B_i| |B_j|)
%
%   The per-recording scalar C is the mean of C_ij over all pairs with
%   |B_i| >= minBursts and |B_j| >= minBursts.
%
%   The companion-of-nearest-neighbour counting convention (each burst in
%   B_i matches its nearest neighbour in B_j and is counted once if that
%   nearest distance is <= windowMs) is the same as used by the tmp/
%   prototype tmp/m/tmp_angle_burstsync.m line 88-97.
%
% INPUTS:
%   burstStartTimes  -  1 x nCh cell of burst start-time vectors in seconds.
%                       Channels with fewer than minBursts bursts are
%                       excluded from the average.
%
% Name-value options (defaults from project_config().burst_sync):
%   'windowMs'   -  coincidence window in ms (default cfg.burst_sync.window_ms = 50;
%                   Chiappalone 2006: burst onsets propagate with < 100 ms jitter)
%   'minBursts'  -  minimum bursts per channel to include a channel in
%                   the pairwise average (default 3)
%
% OUTPUTS:
%   result.C             -  scalar mean pairwise coincidence in [0, 1]
%   result.windowMs      -  echoed window
%   result.minBursts     -  echoed minBursts
%   result.channelsUsed  -  1 x k vector of channel indices passing the
%                            minBursts gate
%   result.nPairs        -  number of pair values averaged
%
% References:
%   Brofiga M, Poli D, Massobrio P, et al. On the functional role of
%     excitatory and inhibitory populations on the overall network
%     dynamics. (MCS 60-ch MEA). 2023. (Burst propagation / initiator
%     cluster analysis.)
%   Chiappalone M, Bove M, Vato A, Tedesco M, Martinoia S. Dissociated
%     cortical networks show spontaneously correlated activity patterns
%     during in vitro development. Brain Res 2006;1093:41-53.
%     (Burst onsets propagate across the array with < 100 ms jitter.)
%
% See also: TRANSFER_ENTROPY_D1, CONNECTIVITY_XCORR.

    %%% ----- parse inputs --------------------------------------------------
    cfg = project_config();
    defWindowMs  = pick_default(cfg, 'burst_sync', 'window_ms', 50);
    defMinBursts = 3;                                 % tmp/ prototype convention

    p = inputParser;
    addRequired(p,  'burstStartTimes', @iscell);
    addParameter(p, 'windowMs',  defWindowMs,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'minBursts', defMinBursts, @(x) isscalar(x) && x >= 1);
    parse(p, burstStartTimes, varargin{:});
    opt = p.Results;

    nCh = numel(burstStartTimes);
    tau = opt.windowMs / 1000;

    %%% ----- channel eligibility ------------------------------------------
    nb = zeros(nCh, 1);
    bst = cell(nCh, 1);
    for k = 1:nCh
        v = burstStartTimes{k};
        if isempty(v)
            bst{k} = zeros(0, 1);
            nb(k) = 0;
        else
            bst{k} = v(:);
            nb(k) = numel(v);
        end
    end
    valid = find(nb >= opt.minBursts);
    if numel(valid) < 2
        result.C            = 0;
        result.windowMs     = opt.windowMs;
        result.minBursts    = opt.minBursts;
        result.channelsUsed = valid(:).';
        result.nPairs       = 0;
        return;
    end

    %%% ----- pairwise coincidence -----------------------------------------
    pairsC = nchoosek(valid, 2);
    vals = nan(size(pairsC, 1), 1);
    for ii = 1:size(pairsC, 1)
        ia = pairsC(ii, 1);
        ib = pairsC(ii, 2);
        a = bst{ia}; b = bst{ib};
        if isempty(a) || isempty(b)
            vals(ii) = 0;
            continue;
        end
        % For each burst in a, find nearest burst in b; count if within tau.
        [~, jj] = min(abs(a - b.'), [], 2);
        d = abs(a - b(jj));
        coin = sum(d <= tau);
        vals(ii) = coin / sqrt(nb(ia) * nb(ib));
    end
    C = mean(vals, 'omitnan');
    if ~isfinite(C); C = 0; end

    %%% ----- pack ---------------------------------------------------------
    result.C            = C;
    result.windowMs     = opt.windowMs;
    result.minBursts    = opt.minBursts;
    result.channelsUsed = valid(:).';
    result.nPairs       = numel(vals);
end

% =========================================================================
function v = pick_default(cfg, section, field, fallback)
    if isfield(cfg, section) && isfield(cfg.(section), field)
        v = cfg.(section).(field);
    else
        v = fallback;
    end
end
