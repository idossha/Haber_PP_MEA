function run_preprocess(varargin)
%RUN_PREPROCESS Batch entrypoint for pipeline/preprocess_and_save.m.
%
%   Called from scripts/run.sh. Accepts the same name-value arguments
%   as pipeline/preprocess_and_save.m ('overwrite', 'study').
%
%   Example:
%     matlab -nodisplay -nosplash -nodesktop -batch \
%         "addpath(genpath('src'));addpath('scripts'); ...
%          run_preprocess('study','all','overwrite',false)"

    fprintf('\n==== run_preprocess ====\n');
    fprintf('cwd: %s\n\n', pwd);
    t0 = tic;
    preprocess_and_save(varargin{:});
    fprintf('\nrun_preprocess done in %.1f s\n', toc(t0));
end
