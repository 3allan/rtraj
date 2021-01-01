function pathToFile = getPathToCache(file)
% 
% Matt Werner (m.werner@vt.edu) - Dec 30, 2020
% 
% Obtain the full path to the indicated file that should be located in the
% cache. If the file should be in the cache but is not yet saved in the
% cache, then the file will be saved into the cache during the current call
% to rtraj.m. Obtaining the full path may sometimes be necessary on older
% versions of Matlab in which paths must be specified in accordance with
% the native operating system (OS).
% 
%    Inputs:
% 
%              file - The indicated file in the cache whose full path is to
%                     be determined in accordance with this computer's
%                     native OS.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%    Outputs:
% 
%        pathToFile - The full path to the desired file in the cache
%                     according to the format used by this computer's
%                     native OS.
% 

% Ensure the specified file is a string/char
checkInput(file)

% Obtain the path that points to where the file is in the cache
pathToMain = pwd;
if (ispc)
    relativePathToFile = strcat('\cache\', file);
elseif (ismac || isunix)
    relativePathToFile = strcat('/cache/', file);
else
    error("Platform not supported.")
end
pathToFile = strcat(pathToMain, relativePathToFile);