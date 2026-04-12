function root = get_project_root()
%GET_PROJECT_ROOT Resolve the absolute path to the Haber_PP_MEA repository root.
%
%   root = GET_PROJECT_ROOT() returns the absolute path of the project root
%   directory by walking up two levels from the location of this file
%   (src/utils/get_project_root.m). This makes every script independent of
%   the user's current working directory or login.
%
% INPUTS:
%   (none)
%
% OUTPUTS:
%   root  -  Character vector with the absolute path to the project root.

thisFile = mfilename('fullpath');
utilsDir = fileparts(thisFile);          % .../src/utils
srcDir   = fileparts(utilsDir);          % .../src
root     = fileparts(srcDir);            % .../Haber_PP_MEA

end
