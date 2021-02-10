function x = cosspace(d1, d2, n)
% 
% Matt Werner (m.werner@vt.edu) - Feb 7, 2021
% 
% Provide cosine spacing for sampling points.
% 
%    Inputs: 
% 
%                d1 - First element of the sample points.
%                     Size: 1-by-1 (scalar)
%                     Units: ?
% 
%                d2 - Final element of the sample points.
%                     Size: 1-by-1 (scalar)
%                     Units: ?
% 
%                 n - Optional(!) Specifies how many sample points to
%                     include.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 

% Enforce that either 2 or 3 inputs are given
narginchk(2, 3)

% Define a linearly-spaced angular coordinate
if (nargin == 2)
    t = linspace(0, pi);
else
    t = linspace(0, pi, n);
end

% Calculate the cosine-spaced sample points
x = d1 + (d2 - d1) * (1 - cos(t)) / 2;