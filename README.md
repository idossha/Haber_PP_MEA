# Haber_PP_MEA

Analysis pipeline for:

**Acute DOI-induced changes in spontaneous cortical network activity: a multi-electrode array study with ketanserin pharmacological verification**

Ido Haber, Ilhan Bok, Benjamin Kutler, Maya Norenberg, Matthew I. Banks, Giulio Tononi, Aviad Hai

*Target journal: Journal of Neural Engineering*

## Overview

This repository contains the MATLAB analysis pipeline for characterising the acute effects of the 5-HT2A agonist DOI on whole-network electrophysiology in dissociated rat cortical cultures recorded on 60-electrode multi-electrode arrays (MEAs). A paired ketanserin (5-HT2A antagonist) control arm establishes receptor dependence.

The pipeline processes raw Tucker-Davis Technologies (TDT) recordings through filtering, spike detection, burst detection, functional connectivity estimation, and graph-theoretic network analysis, producing publication-ready figures and statistical summaries.

## Repository structure

```
├── src/
│   ├── config/
│   │   └── project_config.m        # All parameters, paths, dataset lists
│   ├── pipeline/
│   │   └── preprocess_and_save.m   # Raw TDT → cached spike/burst data
│   ├── analysis/
│   │   ├── connectivity_xcorr.m    # Pairwise cross-correlation adjacency
│   │   ├── network_metrics.m       # Clustering, path length, modularity, σ
│   │   ├── paired_stats.m          # Hierarchical bootstrap + Wilcoxon
│   │   ├── bh_fdr.m               # Benjamini-Hochberg FDR correction
│   │   ├── hedges_g_av.m          # Paired Hedges' g effect size
│   │   └── ...                     # See src/README.md for full listing
│   ├── figures/
│   │   ├── fig_spike_rate.m        # Paired firing-rate plots
│   │   ├── fig_connectivity_summary.m  # Network metric summary panels
│   │   ├── fig_stats_bootstrap.m   # Hierarchical bootstrap CI/p-values
│   │   └── ...                     # 16 figure scripts total
│   └── utils/                      # Shared helpers (filtering, I/O, plotting)
├── scripts/
│   ├── run.sh                      # Batch entrypoint (preprocess → figures)
│   ├── smoke_test.m                # Reproducibility smoke test
│   └── ...
├── authors.json                    # Structured author/affiliation metadata
└── requirements.txt                # MATLAB toolbox dependencies (informational)
```

## Quick start

### Prerequisites

- MATLAB R2020a or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- [TDT MATLAB SDK](https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/) placed at `TDTMatlabSDK/` in the project root

### Usage

```matlab
% Add the source tree to the MATLAB path
addpath(genpath('src'));

% Build the spike/burst cache from raw TDT recordings (run once)
preprocess_and_save();

% Generate figures
fig_spike_rate('doi');
fig_burst_rate('doi');
fig_connectivity_summary('doi');
fig_stats_bootstrap('doi');

% Run the full pipeline via shell
% ./scripts/run.sh
```

All parameters (filter cutoffs, spike threshold, burst ISI criteria, dataset paths) are defined in a single configuration file: [`src/config/project_config.m`](src/config/project_config.m).

See [`src/README.md`](src/README.md) for detailed documentation of every module.

## Key methods

| Step | Method | Citation |
|------|--------|----------|
| Spike detection | 5σ threshold, 1 ms refractory | Schroeter et al. 2015; Brofiga et al. 2023 |
| Burst detection | ISI state machine (100/200 ms, min 5 spikes) | Chiappalone et al. 2006; Wagenaar et al. 2006 |
| Functional connectivity | Pearson cross-correlation, ±100 ms lag | Garofalo et al. 2009 |
| Network topology | Clustering, path length, modularity Q, σ | Rubinov & Sporns 2010 |
| Statistical inference | Two-level hierarchical bootstrap + BH-FDR | Saravanan et al. 2020 |

## Data availability

Raw TDT recording blocks are available from the corresponding author (A. Hai, ahai@wisc.edu) on reasonable request.

## License

MIT

## Citation

If you use this code, please cite:

> Haber I, Bok I, Kutler B, Norenberg M, Banks MI, Tononi G, Hai A. Acute DOI-induced changes in spontaneous cortical network activity: a multi-electrode array study with ketanserin pharmacological verification. *Journal of Neural Engineering* (in preparation). 2026.
