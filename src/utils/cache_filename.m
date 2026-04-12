function name = cache_filename(datasetName)
%CACHE_FILENAME Build a filesystem-safe cache filename from a dataset name.
%
%   name = CACHE_FILENAME(datasetName) returns a string of the form
%   "cache_<safeName>.mat" where unsafe characters in datasetName ('-' and
%   '#') have been replaced with underscores. This is the canonical naming
%   used by the preprocessing pipeline and every figure script.
%
% INPUTS:
%   datasetName  -  Dataset folder name (e.g. 'IdoDOI-230914-142502_#1').
%
% OUTPUTS:
%   name  -  Cache filename (e.g. 'cache_IdoDOI_230914_142502__1.mat').

    safe = strrep(strrep(datasetName, '-', '_'), '#', '_');
    name = sprintf('cache_%s.mat', safe);
end
