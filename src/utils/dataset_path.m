function dpath = dataset_path(datasetName, cfg)
%DATASET_PATH Resolve a dataset name to its absolute folder on disk.
%
%   dpath = DATASET_PATH(datasetName, cfg) returns the absolute path of the
%   TDT block directory for DATASETNAME by checking the configured DOI and
%   ketanserin data folders. The function errors if no matching folder is
%   found in either location.
%
% INPUTS:
%   datasetName  -  Folder name (e.g. 'IdoControl-230914-130200_#1').
%   cfg          -  Project config struct from project_config().
%
% OUTPUTS:
%   dpath  -  Absolute path to the TDT block directory.

    candidates = { ...
        fullfile(cfg.paths.data_doi, datasetName), ...
        fullfile(cfg.paths.data_ket, datasetName)};
    for k = 1:numel(candidates)
        if exist(candidates{k}, 'dir')
            dpath = candidates{k};
            return;
        end
    end
    error('dataset_path:NotFound', ...
        ['Dataset folder "%s" not found in either data root:\n', ...
         '  %s\n  %s'], datasetName, candidates{1}, candidates{2});
end
