function gen_standalone_adj()
%GEN_STANDALONE_ADJ Export adjacency matrices as CSVs for Python plotting.
%
%   GEN_STANDALONE_ADJ() computes cross-correlation adjacency matrices for
%   DOI and ketanserin exemplar pairs and writes them as CSV files to
%   figures/panels/standalone_adj/.
%
%   Then run the Python plotter:
%       python3 scripts/plot_standalone_adj.py
%
%   Exemplar pair indices match the main-text figure choices:
%       DOI connectivity: pair 2
%       Ket connectivity: pair 2
%
% See also: FIG_CONNECTIVITY_EXEMPLAR, CONNECTIVITY_XCORR.

    cfg = project_config();

    studies = struct( ...
        'name',      {'doi',  'ket'}, ...
        'pairIndex', {2,      2});

    for s = 1:numel(studies)
        study = studies(s).name;
        pairIdx = studies(s).pairIndex;

        panelDir = output_path(cfg, study, 'connectivity', 'panels');
        if ~exist(panelDir, 'dir'); mkdir(panelDir); end

        fprintf('\n=== %s (pair %d) ===\n', upper(study), pairIdx);

        [pairs, labels] = get_pairs_and_labels(cfg, study);
        pair = pairs(pairIdx);

        fprintf('  baseline : %s\n', pair.baseline);
        fprintf('  treatment: %s\n', pair.treatment);

        channels = cfg.channels.default;
        [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, channels, cfg);

        xcorrArgs = { ...
            'binMs',         cfg.connectivity.bin_ms, ...
            'maxLagMs',      cfg.connectivity.max_lag_ms, ...
            'normalization', cfg.connectivity.normalization};
        bRes = connectivity_xcorr(bSpikes, 'durationSec', bMeta.durationSec, xcorrArgs{:});
        tRes = connectivity_xcorr(tSpikes, 'durationSec', tMeta.durationSec, xcorrArgs{:});
        delta = tRes.adjacency - bRes.adjacency;

        % Save as CSV
        writematrix(bRes.adjacency, fullfile(panelDir, [study '_adj_baseline.csv']));
        writematrix(tRes.adjacency, fullfile(panelDir, [study '_adj_treatment.csv']));
        writematrix(delta,          fullfile(panelDir, [study '_adj_delta.csv']));

        fprintf('  Wrote 3 CSVs to %s\n', panelDir);
    end

    fprintf('\n=== Done. Now run: python3 scripts/plot_standalone_adj.py ===\n');
end
