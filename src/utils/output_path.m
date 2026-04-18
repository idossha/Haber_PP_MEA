function dir = output_path(cfg, study, category, subdir)
%OUTPUT_PATH Resolve output directory for a given study and figure category.
%
%   dir = OUTPUT_PATH(cfg, 'doi', 'rates', 'panels')
%       → fullfile(cfg.paths.output, 'fig2', 'panels')
%
%   dir = OUTPUT_PATH(cfg, 'ket', 'connectivity', 'stats')
%       → fullfile(cfg.paths.output, 'fig5', 'stats')
%
%   dir = OUTPUT_PATH(cfg, '', 'si', '')
%       → fullfile(cfg.paths.output, 'SI')
%
%   dir = OUTPUT_PATH(cfg, '', 'shared', '')
%       → fullfile(cfg.paths.output, 'shared')
%
% INPUTS:
%   cfg      - project_config() struct (needs cfg.paths.output)
%   study    - 'doi' or 'ket' (ignored for 'si' and 'shared')
%   category - 'rates', 'connectivity', 'si', or 'shared'
%   subdir   - 'panels', 'stats', 'exemplars', 'composite', or ''
%
% OUTPUT:
%   dir - absolute path to the resolved directory

    figMap = struct( ...
        'doi_rates',        'fig2', ...
        'doi_connectivity', 'fig3', ...
        'ket_rates',        'fig4', ...
        'ket_connectivity', 'fig5');

    switch lower(category)
        case 'si'
            dir = fullfile(cfg.paths.output, 'SI');
        case 'shared'
            dir = fullfile(cfg.paths.output, 'shared');
        otherwise
            key = sprintf('%s_%s', lower(study), lower(category));
            if ~isfield(figMap, key)
                error('output_path:BadKey', ...
                    'Unknown study/category combination: %s/%s', study, category);
            end
            figName = figMap.(key);
            if isempty(subdir)
                dir = fullfile(cfg.paths.output, figName);
            else
                dir = fullfile(cfg.paths.output, figName, subdir);
            end
    end
end
