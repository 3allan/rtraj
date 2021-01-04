function remainingLines = getNumberOfLinesLeft(fid, lineNumber)
% 
% Matt Werner (m.werner@vt.edu) - Jan 4, 2020
% 
% Count the remaining amount of lines in a file from the current pointer
% position by iterating through the remainder of the file with file ID fid
% and tracking how many newline characters (\n) there are. Also, reset the
% pointer back to where it was on the function call.
% 
% 
%    Inputs:
% 
%               fid - File ID of the current file being read.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%        lineNumber - Optional(!) Current line number that the pointer is
%                     on in the file. Leaving this quantity unspecified
%                     assumes that the pointer's current position is at the
%                     beginning (where the pointer is placed upon opening
%                     the file).
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%    Outputs:
% 
%    remainingLines - Number of remaining lines from the current position
%                     of the pointer.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 

% Count the number of lines remaining
remainingLines = getNumberOfLines(fid);

% Restore the pointer to its original location
repositionPointer(fid, lineNumber)