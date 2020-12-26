function checkInput(varargin)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Check if the input s is of type string (or character)).
% 
%    Inputs:
% 
%                 s - Object to check type. Must be either of type (class)
%                     <string> or <char> so that no error is produced.
%                     Size: 1-by-1 (object)
%                     Units: N/A
% 
%    Outputs:
% 
%                   -
% 

numelargin = numel(varargin);
for ii = 1:numelargin
    s = varargin{ii};
    if (~(ischar(s) || isstring(s)))
        error("Please provide a model as a string ("""") or character ('').")
    end
end