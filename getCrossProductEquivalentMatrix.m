function xProdEquivMatrix = getCrossProductEquivalentMatrix(v)
% 
% Matt Werner (m.werner@vt.edu) - Dec 23, 2020
% 
% Form the cross-product equivalent matrix from the 3-vector v such that,
% given another 3-vector, the product cross(v, u) is equivalent to
% performing matrix multiplication. Explicitly, the cross-product
% equivalent matrix is given as
%            _                    _
%           |  0       -v       v  |
%           |            3       2 |
%           |                      |
%    [v*] = |                      | ,
%           |  v        0      -v  |
%           |   3                1 |
%           |                      |
%           |                      |
%           | -v        v        0 |
%           |_  2        1        _|
% 
% where [v*] is the cross-product equivalent matrix such that 
% cross(v, u) = [v*]*u.
% 
%    Inputs:
% 
%                 v - Column vector containing 3 elements to be placed into
%                     the cross-product equivalent matrix form according to
%                     its shown definition.
%                     Size: 3-by-1 (vector)
%                     Units: ?
% 
%    Outputs:
% 
%  xProdEquivMatrix - Cross-product equivalent matrix of the given 3-vector.
%                     Size: 3-by-3 (matrix)
%                     Units: ?
% 

% Check that the given vector exists in 3D
if (~all(size(v(:)) == [3, 1]))
    error("Vector must exist in 3D.")
end

% Explicitly designate components
v1 = v(1);
v2 = v(2);
v3 = v(3);

% Create the cross-product equivalent matrix
xProdEquivMatrix = [0, -v3, v2; v3, 0, -v1; -v2, v1, 0];