function paths = export_figure_stats(stats, basePath)
%EXPORT_FIGURE_STATS Write a figure's numeric summary as CSV + JSON sidecars.
%
%   paths = EXPORT_FIGURE_STATS(stats, basePath) serialises the flat
%   struct STATS (scalar or 1x1) to two companion files next to a
%   figure:
%       <basePath>.csv   - single-row key/value table
%       <basePath>.json  - human-readable JSON
%
%   The intent is that every figure script calls this once with the
%   numbers it computed, so the paper writer has a single location
%   (figures_out/*.csv or *.json) from which to pull every in-text
%   number without opening MATLAB again.
%
% INPUTS:
%   stats     -  Struct of scalar / short-vector numeric fields (strings
%                and char vectors are also OK). Vectors of length >= 2
%                are expanded into `<field>_1`, `<field>_2`, ...
%   basePath  -  Output path without extension. Parent directory is
%                created if it does not exist.
%
% OUTPUTS:
%   paths  -  Struct with fields .csv and .json (absolute paths written).
%
% See also: PAIRED_STATS, FIG_SPIKE_RATE, FIG_BURST_RATE,
%           FIG_CONNECTIVITY_SUMMARY.

    if ~isstruct(stats)
        error('export_figure_stats:Args', 'stats must be a struct.');
    end
    parentDir = fileparts(basePath);
    if ~isempty(parentDir) && ~exist(parentDir, 'dir')
        mkdir(parentDir);
    end

    paths.csv  = [basePath '.csv'];
    paths.json = [basePath '.json'];

    % --- Flatten short vectors into scalar fields ----------------------
    flat = struct();
    fns  = fieldnames(stats);
    for k = 1:numel(fns)
        name = fns{k};
        val  = stats.(name);
        if isnumeric(val) && ~isscalar(val) && numel(val) <= 8
            for i = 1:numel(val)
                flat.(sprintf('%s_%d', name, i)) = val(i);
            end
        elseif ischar(val) || isstring(val)
            flat.(name) = char(val);
        elseif isnumeric(val) && isscalar(val)
            flat.(name) = double(val);
        elseif islogical(val) && isscalar(val)
            flat.(name) = double(val);
        elseif isnumeric(val)
            % Larger vectors / matrices: store as JSON only.
            flat.(name) = val;
        else
            % Unknown types become their text representation.
            try
                flat.(name) = char(string(val));
            catch
                flat.(name) = '';
            end
        end
    end

    % --- CSV (flat scalars only) ---------------------------------------
    fid = fopen(paths.csv, 'w');
    if fid == -1
        error('export_figure_stats:WriteFailed', ...
            'Could not open %s for writing.', paths.csv);
    end
    cleanup = onCleanup(@() fclose(fid));
    csvFields = fieldnames(flat);
    scalarMask = false(numel(csvFields), 1);
    for k = 1:numel(csvFields)
        v = flat.(csvFields{k});
        scalarMask(k) = isnumeric(v) && isscalar(v) ...
                     || islogical(v) && isscalar(v) ...
                     || ischar(v);
    end
    keep = csvFields(scalarMask);
    fprintf(fid, '%s\n', strjoin(keep, ','));
    vals = cell(1, numel(keep));
    for k = 1:numel(keep)
        v = flat.(keep{k});
        if ischar(v)
            vals{k} = csv_escape(v);
        elseif isnan(v)
            vals{k} = 'NaN';
        elseif isinf(v)
            if v > 0; vals{k} = 'Inf'; else; vals{k} = '-Inf'; end
        else
            vals{k} = sprintf('%.10g', v);
        end
    end
    fprintf(fid, '%s\n', strjoin(vals, ','));
    clear cleanup;

    % --- JSON (full struct, matrices included) -------------------------
    try
        txt = jsonencode(flat, 'PrettyPrint', true);
    catch
        % Older MATLAB without 'PrettyPrint'.
        txt = jsonencode(flat);
    end
    fid = fopen(paths.json, 'w');
    if fid == -1
        error('export_figure_stats:WriteFailed', ...
            'Could not open %s for writing.', paths.json);
    end
    fwrite(fid, txt, 'char');
    fclose(fid);
end

% -------------------------------------------------------------------------
function s = csv_escape(txt)
    if any(txt == ',' | txt == '"' | txt == newline)
        s = ['"' strrep(txt, '"', '""') '"'];
    else
        s = txt;
    end
end
