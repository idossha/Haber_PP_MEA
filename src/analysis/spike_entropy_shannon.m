function result = spike_entropy_shannon(spikeCells, durationSec, varargin)
%SPIKE_ENTROPY_SHANNON Per-channel Shannon entropy of binned spike counts.
%
%   result = SPIKE_ENTROPY_SHANNON(spikeCells, durationSec) computes the
%   per-channel Shannon entropy H = -sum(p log2 p) of binned spike-count
%   histograms, with counts clipped at cfg.entropy.clip_counts to yield a
%   stable discrete alphabet (Varley 2024 §2.2 convention).
%
%   Note: entropy production (time-reversibility) is NOT reported here.
%   The tmp/NEW_ANGLES_REPORT.md §2 caveats section documents that EP is
%   too noisy at our recording length (~10 min) and has been dropped
%   from the paper. Only the marginal Shannon H survives as a per-well
%   scalar for paired_stats.
%
% INPUTS:
%   spikeCells    -  1 x nCh cell of column-vector spike times (seconds)
%   durationSec   -  scalar recording duration (seconds)
%
% Name-value options (defaults from project_config().entropy):
%   'binMs'       -  bin width in ms (default cfg.entropy.bin_ms = 50;
%                     Varley 2024 §2.2 50 ms convention)
%   'clipCounts'  -  clip value for discrete alphabet (default
%                     cfg.entropy.clip_counts = 3; Varley 2024 §2.2 —
%                     yields alphabet {0,1,2,3} with H_max = 2 bits)
%
% OUTPUTS:
%   result.H           -  nCh x 1 per-channel Shannon entropy in bits;
%                          NaN for silent or empty channels
%   result.medianH     -  scalar, median of H across channels (NaN-safe)
%   result.binMs       -  echoed
%   result.clipCounts  -  echoed
%
% References:
%   Varley TF, Pope M, Grazzi LF, Pittman-Polletta B, Sporns O.
%     Information processing dynamics in neural networks of macroscopic
%     brain organoids. J Neural Eng 2024. (§2.2 — 50 ms bin, clip-3
%     Shannon H convention.)
%
% See also: TRANSFER_ENTROPY_D1, BURST_ONSET_COINCIDENCE.

    %%% ----- parse inputs --------------------------------------------------
    cfg = project_config();
    defBinMs      = pick_default(cfg, 'entropy', 'bin_ms',      50);
    defClipCounts = pick_default(cfg, 'entropy', 'clip_counts', 3);

    p = inputParser;
    addRequired(p,  'spikeCells',  @iscell);
    addRequired(p,  'durationSec', @(x) isscalar(x) && isfinite(x) && x > 0);
    addParameter(p, 'binMs',      defBinMs,      @(x) isscalar(x) && x > 0);
    addParameter(p, 'clipCounts', defClipCounts, @(x) isscalar(x) && x >= 1);
    parse(p, spikeCells, durationSec, varargin{:});
    opt = p.Results;

    nCh = numel(spikeCells);
    H = nan(nCh, 1);

    %%% ----- bin edges -----------------------------------------------------
    edges = 0:(opt.binMs / 1000):durationSec;
    nBins = numel(edges) - 1;
    if nBins < 10
        result.H          = H;
        result.medianH    = NaN;
        result.binMs      = opt.binMs;
        result.clipCounts = opt.clipCounts;
        return;
    end

    %%% ----- per-channel Shannon entropy ----------------------------------
    % Varley 2024 §2.2: discrete alphabet {0, 1, ..., clipCounts}; counts
    % above clipCounts are clipped to clipCounts. H = -sum p log2 p over
    % the marginal distribution of bin values.
    for c = 1:nCh
        st = spikeCells{c};
        if isempty(st)
            H(c) = 0;
            continue;
        end
        x = histcounts(st, edges);
        x(x > opt.clipCounts) = opt.clipCounts;
        vals = unique(x);
        q = histcounts(x, [vals, max(vals) + 1]) / nBins;
        q = q(q > 0);
        if isempty(q)
            H(c) = 0;
        else
            H(c) = -sum(q .* log2(q));
        end
    end

    %%% ----- pack ---------------------------------------------------------
    result.H          = H;
    result.medianH    = median(H, 'omitnan');
    result.binMs      = opt.binMs;
    result.clipCounts = opt.clipCounts;
end

% =========================================================================
function v = pick_default(cfg, section, field, fallback)
    if isfield(cfg, section) && isfield(cfg.(section), field)
        v = cfg.(section).(field);
    else
        v = fallback;
    end
end
