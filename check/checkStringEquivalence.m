function checkStringEquivalence(strArray, flag)
% 
% Matt Werner (m.werner@vt.edu) - Jan 7, 2021
% 
% Check if all strings contained within the provided array of strings are
% equivalent to one another. The only exception to absolute equivalence
% between all strings in the array is the empty string, if desired.
% 
%    Inputs:
% 
%          strArray - Array of strings whose elements are to be compared
%                     against one another for string equivalence.
%                     Size: n-by-m (string matrix)
%                     Units: N/A
% 
%              flag - Optional(!) Indication (true/false) as to whether
%                     empty strings should be included in the determination
%                     of string array equivalence. That is, empty strings
%                     are included in comparisons against other strings in
%                     the array if true and are otherwise not included if
%                     false. The DEFAULT behavior is TRUE (to include empty
%                     strings in comparisons).
%                     Size: 1-by-1 (boolean)
%                     Units: N/A
% 
%    Outputs:
% 
%                   -
% 

% Enforce that the flag is a boolean
if (~islogical(flag))
    error("Flag must be of type boolean.")
end

% Get size of the string array
numelStrArray = numel(strArray);

% Cannot fail the check if only 1 string is supplied
if (numelStrArray == 1)
    return
end

% Define default behavior
flag_includeEmptyStrings = true;
if (nargin == 2 || flag == false)
    flag_includeEmptyStrings = false;
end
    

% Perform the checks element-by-element against the first element,
% skipping self-equivalence
s = strArray(1);
for ii = 2:numelStrArray
    % Skip if empty strings shouldn't be included in the check
    if (~flag_includeEmptyStrings), continue, end
    % Compare the (non-empty) strings
    equivalentStrings = strcmp(s, strArray(ii));
    if (~equivalentStrings)
        % Get the matrix element that failed
        strArrayNumberOfRows = size(strArray, 1);
        col = ceil(ii/strArrayNumberOfRows);
        row = ii - (strArrayNumberOfRows * (col - 1));
        % Send error message
        error("Mismatched strings (%1.0f, %1.0f).", row, col)
    end
end