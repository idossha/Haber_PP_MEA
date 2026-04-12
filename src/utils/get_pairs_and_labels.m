function [pairs, labels] = get_pairs_and_labels(cfg, study)
%GET_PAIRS_AND_LABELS Build the dataset pair list and axis labels for a study.
%
%   [pairs, labels] = GET_PAIRS_AND_LABELS(cfg, study) returns the paired
%   baseline/treatment dataset list (from pair_datasets()) and a struct
%   with .baseline / .treatment label strings for the requested study.
%
% INPUTS:
%   cfg    -  Project config struct from project_config().
%   study  -  'doi' or 'ket' (case-insensitive).
%
% OUTPUTS:
%   pairs   -  Struct array from pair_datasets().
%   labels  -  Struct with .baseline / .treatment label strings.

    study = lower(study);
    switch study
        case 'doi'
            pairs  = pair_datasets(cfg.datasets.doi.baseline, cfg.datasets.doi.treatment);
            labels = cfg.labels.doi;
        case 'ket'
            pairs  = pair_datasets(cfg.datasets.ket.baseline, cfg.datasets.ket.treatment);
            labels = cfg.labels.ket;
        otherwise
            error('get_pairs_and_labels:UnknownStudy', ...
                'study must be one of {''doi'', ''ket''}, got "%s".', study);
    end
end
