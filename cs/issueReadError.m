function issueReadError(msg, file, lineNumber)
% 
% Matt Werner (m.werner@vt.edu) - Jan 1, 2021
% 
% Issue an error when reading input files into rtraj.
% 
%    Inputs:
%               msg - Error message.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%              file - Path to the file where in the errors occurs.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%        lineNumber - Line number on which the error occurs.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
%                   
 
errorNameAndLine = " (%s: line %1.0f)";
error(msg + errorNameAndLine, file, lineNumber)