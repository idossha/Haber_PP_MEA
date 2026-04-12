function channels = select_representative_channels(baselineRates, treatmentRates, nWanted, varargin)
%SELECT_REPRESENTATIVE_CHANNELS Pick SNR + story-consistent channels for a trace plot.
%
%   channels = SELECT_REPRESENTATIVE_CHANNELS(baselineRates, treatmentRates, nWanted)
%   chooses NWANTED channel indices optimised for (a) good baseline
%   signal-to-noise (channels with detectable spiking, not flat or
%   saturated), and (b) the "study story" (channels where the treatment
%   changes rate in the direction of the group-level effect). The result
%   is a sorted 1 x NWANTED row vector of channel indices.
%
%   Strategy:
%     1. Eligibility filter:
%           baselineRates >= rateFloor        (not dead on baseline)
%           baselineRates <= rateCeilingPct   (not rail-saturated)
%           treatmentRates >= 0.05 (alive on treatment)
%     2. Story score (per channel):
%           scoreStory(c) = (treatment - baseline) / max(1, baseline)
%        — signed direction. The group-level sign (positive for DOI,
%        negative for ketanserin) is inferred from the median of
%        scoreStory across eligible channels, then the score is flipped
%        so "high score" = "follows the story".
%     3. Diversity: from the top-ranked story channels, pick NWANTED
%        that evenly span the baseline spike-rate percentiles (so the
%        resulting panel shows a range of rate levels, not just the
%        most-active channels).
%
%   SELECT_REPRESENTATIVE_CHANNELS(..., 'Name', value, ...) options:
%     'rateFloor'            minimum baseline spikes/min (default 5)
%     'rateCeilingPercentile' upper bound percentile (default 98)
%     'topStoryFraction'     fraction of eligible channels kept after
%                            ranking by the story score (default 0.7 =
%                            keep the top 70 %)
%
% INPUTS:
%   baselineRates   -  (nCh x 1) per-channel baseline rate (spikes/min)
%   treatmentRates  -  (nCh x 1) per-channel treatment rate (spikes/min)
%   nWanted         -  Number of channels to return (e.g. 10)
%
% OUTPUTS:
%   channels  -  Sorted row vector of channel indices.

    p = inputParser;
    addRequired(p,  'baselineRates');
    addRequired(p,  'treatmentRates');
    addRequired(p,  'nWanted', @(x) isscalar(x) && x > 0);
    addParameter(p, 'rateFloor',             5.0);
    addParameter(p, 'rateCeilingPercentile', 98);
    addParameter(p, 'topStoryFraction',      0.7);
    parse(p, baselineRates, treatmentRates, nWanted, varargin{:});
    opt = p.Results;

    b = baselineRates(:);
    t = treatmentRates(:);
    nCh = numel(b);
    if numel(t) ~= nCh
        error('select_representative_channels:Args', ...
            'baselineRates and treatmentRates must be the same length.');
    end

    % --- Eligibility: alive, not saturated -----------------------------
    alive = b >= opt.rateFloor & t >= 0.05;
    if any(alive)
        ceil_r = prctile(b(alive), opt.rateCeilingPercentile);
        keep   = alive & b <= ceil_r;
    else
        keep   = false(nCh, 1);
    end
    eligible = find(keep);

    if numel(eligible) <= nWanted
        % Not enough candidates — return what we have.
        channels = sort(eligible(:))';
        return;
    end

    % --- Story score: signed relative change ---------------------------
    bE = b(eligible);
    tE = t(eligible);
    scoreStory = (tE - bE) ./ max(1, bE);

    % Infer group-level direction (positive for activation, negative for
    % suppression), then flip so higher = more on-story.
    groupSign = sign(median(scoreStory));
    if groupSign == 0; groupSign = 1; end
    scoreOnStory = scoreStory * groupSign;

    % --- Keep top X % by story score -----------------------------------
    nKeepTop = max(nWanted, round(numel(eligible) * opt.topStoryFraction));
    [~, orderStory] = sort(scoreOnStory, 'descend');
    topIdxInElig = orderStory(1:min(nKeepTop, numel(eligible)));
    topChanIdx   = eligible(topIdxInElig);
    topBaseRate  = b(topChanIdx);

    % --- Diversity: pick nWanted spanning baseline rate percentiles ----
    [~, orderRate] = sort(topBaseRate, 'ascend');
    rankPicks = round(linspace(1, numel(orderRate), nWanted));
    rankPicks = unique(rankPicks, 'stable');
    while numel(rankPicks) < nWanted
        extra = randi(numel(orderRate));
        if ~ismember(extra, rankPicks)
            rankPicks(end+1) = extra; %#ok<AGROW>
        end
    end
    chosen = topChanIdx(orderRate(rankPicks));
    channels = sort(chosen(:))';
end
