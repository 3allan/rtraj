function restorePointer(fid, requestedLine)
% 
% Matt Werner (m.werner@vt.edu) - Jan 4, 2021
% 
% Reposition the pointer reading through a file to a requested line number
% by counting line numbers from the beginning and iterating through until
% the current line number matches the requested line number.
% 
%    Inputs:
%               fid - File ID of the current file whose pointer is being
%                     repositioned to the requested line number.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%     requestedLine - Optional(!) Requested line number at which the
%                     pointer reading through the file be placed.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 

% Check that the current line number is positive
if (lineNumber <= 0)
    error('Requested line number must be positive.')
end

% Reset the pointer to the beginning
frewind(fid);

% Behavior defined to be finished if current line number is unspecified
if (nargin == 1 || lineNumber == 1), return, end

% Reset the pointer to its original position
currentLine = 1;
while (currentLine <= requestedLine)
    fgetl(fid);
    currentLine = currentLine + 1;
end