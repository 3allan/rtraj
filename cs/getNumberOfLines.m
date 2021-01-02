function n = getNumberOfLines(fid)
% 
% Matt Werner (m.werner@vt.edu) - Jan 2, 2021
% 
% Get the number of (remaining) lines in a file by passing the file ID
% (fid) and counting the number of newline characters (\n) from the current
% pointer position.
% 
%    Inputs:
% 
%               fid - File identifier that specifies the location in the
%                     file being read.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%    Outputs:
%                 n - Number of (remaining) lines in the file. The total
%                     number of lines in the file is obtained (ie n = N)
%                     when the file ID is passed directly after opening the
%                     file up.
% 

% Begin the counter
n = 0;

% Count the number of (remaining) lines the file
while (~feof(fid))
    fgetl(fid);
    n = n + 1;
end