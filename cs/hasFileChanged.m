function profilesToLoad = hasFileChanged(file, modified, ...
                               previous_file, previous_modified)
% 
% Matt Werner (m.werner@vt.edu) - Dec 31, 2020
% 
% Determine which, if any, files have changed either by name or contents.
% The files to be compared here are likely to be ones with very similar
% contents to one another, but differ in only a few lines or by name. A
% file changing its name without changing its contents may be seemingly
% inconsequential, but is counted as a modification to the file in this
% context to maximize confidence that a mistake from ambiguity is avoided
% in reading the file's contents.
% 
%    Inputs:
% 
%              file - The name(s) of the current file(s) holding the most
%                     up-to-date contents.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 
%          modified - The date and times at which the current file(s)
%                     was/were last modified.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 
%     previous_file - The name(s) of the previous file(s) holding the
%                     previously most up-to-date contents.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 
% previous_modified - The date and times at which the previous file(s)
%                     was/were last modified.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 
%    Outputs:
% 
%    profilesToLoad - Indication of which file(s) have changed.
%                     Size: n-by-1 (logical array)
%                     Units: N/A
% 

profilesToLoad = find(~strcmp(file, previous_file) + ...
                      ~strcmp(modified, previous_modified));