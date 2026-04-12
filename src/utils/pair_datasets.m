function pairs = pair_datasets(baselineList, treatmentList)
%PAIR_DATASETS Match baseline and treatment dataset names by "#N" suffix.
%
%   pairs = PAIR_DATASETS(baselineList, treatmentList) returns a struct
%   array, one entry per pair, with fields:
%       .id         - the integer N parsed from the trailing "#N" suffix
%       .baseline   - baseline dataset name
%       .treatment  - matching treatment dataset name
%
%   The two input cell arrays must contain dataset names ending in
%   "_#<number>" (e.g. 'IdoControl-230914-130200_#1'). Pairs are returned in
%   ascending order of N. The function errors if either list contains a
%   name without a "#N" suffix, or if any baseline N has no matching
%   treatment N (or vice versa).
%
% INPUTS:
%   baselineList   -  Cell array of baseline dataset folder names.
%   treatmentList  -  Cell array of treatment dataset folder names.
%
% OUTPUTS:
%   pairs  -  Struct array with .id, .baseline, .treatment.

    if ~iscell(baselineList) || ~iscell(treatmentList)
        error('pair_datasets:Args', ...
            'Both inputs must be cell arrays of dataset names.');
    end

    bIds = parseIds(baselineList,  'baselineList');
    tIds = parseIds(treatmentList, 'treatmentList');

    [commonIds, ibCommon, itCommon] = intersect(bIds, tIds);
    if numel(commonIds) ~= numel(bIds) || numel(commonIds) ~= numel(tIds)
        missingFromT = setdiff(bIds, tIds);
        missingFromB = setdiff(tIds, bIds);
        msg = '';
        if ~isempty(missingFromT)
            msg = sprintf('%sBaseline #%s have no treatment match. ', ...
                msg, sprintfList(missingFromT));
        end
        if ~isempty(missingFromB)
            msg = sprintf('%sTreatment #%s have no baseline match. ', ...
                msg, sprintfList(missingFromB));
        end
        error('pair_datasets:Unmatched', '%s', msg);
    end

    [commonIds, sortIdx] = sort(commonIds);
    ibCommon = ibCommon(sortIdx);
    itCommon = itCommon(sortIdx);

    nPairs = numel(commonIds);
    pairs = repmat(struct('id', [], 'baseline', '', 'treatment', ''), 1, nPairs);
    for k = 1:nPairs
        pairs(k).id        = commonIds(k);
        pairs(k).baseline  = baselineList{ibCommon(k)};
        pairs(k).treatment = treatmentList{itCommon(k)};
    end
end

% -------------------------------------------------------------------------
function ids = parseIds(list, label)
    ids = nan(1, numel(list));
    for k = 1:numel(list)
        tok = regexp(list{k}, '_#(\d+)$', 'tokens', 'once');
        if isempty(tok)
            error('pair_datasets:NoSuffix', ...
                '%s entry "%s" does not end in _#<integer>.', ...
                label, list{k});
        end
        ids(k) = str2double(tok{1});
    end
end

% -------------------------------------------------------------------------
function s = sprintfList(ids)
    s = strjoin(arrayfun(@(n) sprintf('#%d', n), ids, 'UniformOutput', false), ', ');
end
