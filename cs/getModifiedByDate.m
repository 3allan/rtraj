function modifiedByDate = getModifiedByDate(filePath)
% 
% Matt Werner (m.werner@vt.edu) - Dec 30, 2020
% 
% Obtain the modified-by date(s) for the file path(s) indicated by the
% collection of strings contained in filePath. 
% 
%    Inputs:
% 
%          filePath - Path to the file whose modified-by date is sought.
%                     Multiple file paths may be passed in as a vector.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 
%    Outputs:
% 
%    modifiedByDate - Date and time at which the indicated file was last
%                     written to the disk.
%                     Size: n-by-1 (string vector)
%                     Units: N/A
% 

% No checks - Let Matlab errors handle the case that a file name is
% misspelled or otherwise missing

% Allocate memory
modifiedByDate = repmat("", size(filePath));

% Iterate through each file in the indicated path
for ii = 1:numel(filePath)
    filePathii = filePath(ii);
    if (strcmp(filePathii, "")), continue, end
    modifiedByDate(ii) = dir(filePathii).date;
end