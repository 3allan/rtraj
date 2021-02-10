function fq = fastinterp2(x, y, f, xq, yq)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Quickly interpolate on a 2D grid by means of bilinear interpolation. The
% method finds the nearest 4 points to (xq, yq) and performs linear
% interpolation in both the x and y directions according to the formula
%                               [f(x1, y1), f(x1, y2);   [y2 - yq,
%                [x2-xq, xq-x1]  f(x2, y1), f(x2, y2) ]    y - y1 ]
%    f(xq, yq) ~ --------------------------------------------------.
%                                (x2 - x1) (y2 - y1)
% It is assumed that the specified (x, y) are strictly monotone increasing.
% The spacing between each value in both x and y need not be constant since
% only the most local region possible is considered. Special cases like
% repeating values in x or y, boundary points, points outside of the
% domain, etc. are not considered or treated. Thus, it is assumed that
% neither x nor y contain repeating values and that the point (xq, yq) is
% contained inside of the domain (not on or outside of the boundary).
% Note that bilinear interpolation actually results in evaluating a
% quadratic interpolant.
% 
%    Inputs:
% 
%                 x - Determines the x values at which the function f was
%                     evaluated.
%                     Size: 1-by-m (vector)
%                     Units: ?
% 
%                 y - Determines the y values at which the function f was
%                     evaluated.
%                     Size: n-by-1 (vector)
%                     Units: ?
% 
%                 f - The function that was evaluated at (x, y) and is to
%                     be interpolated at the query points (xq, yq).
%                     Size: n-by-m (matrix)
%                     Units: ?
% 
%                xq - Query point in x to evaluate the function f at by
%                     bilinear interpolation.
%                     Size: 1-by-1 (scalar)
%                     Units: ? (same as x)
% 
%                yq - Query point in y to evaluate the function f at by
%                     bilinear interpolation.
%                     Size: 1-by-1 (scalar)
%                     Units: ? (same as y)
% 
%    Outputs:
% 
%                fq - Approximate value of f(xq, yq) calculated by bilinear
%                     interpolation.
%                     Size: 1-by-1 (scalar)
%                     Units: ? (same as f)
% 

% Perform only the most time-efficient and necessary checks
n = max(size(y)); % Amount of rows in y
m = max(size(x)); % Amount of column in x
% Check x and y are at least of dimension 2 to perform the interpolation
if (any(n < 2 | m < 2))
    error("Input elements x and y must be at least of dimension 2.")
end
% Check if sizes of the pair (x, y) and f are compatible
if (~all([n, m] == size(f)))
    error("Sizes of meshgrid(x, y) and f must be the same.")
end
% Check if the input is monotone increasing by checking only the following
% element and none of the others.
if (x(2) == x(1) || y(2) == y(1))
    error("Input elements x and y must be strictly increasing.")
end
% Try to fix if either x or y is monotone decreasing
if (x(2) < x(1))
    % Flip the x direction and the columns in f
    x = flip(x);
    f = flip(f, 2);
end
if (y(2) < y(1))
    % Flip the y direction and the rows in f
    y = flip(y);
    f = flip(f, 1);
end

% Find the closest 4 points surrounding (Xq, Yq) to perform bilinear
% interpolation
L = find(x <= xq, 1, 'last'); % Left index (x direction)
R = find(x >= xq, 1, 'first'); % Right index (x direction)
T = find(y >= yq, 1, 'first'); % Top index (y direction)
B = find(y <= yq, 1, 'last'); % Bottom index (y direction)

% Check if there were any direct matches so that no interpolation is
% necessary
if (L == R && T == B)
    fq = f(B, L); % Arbitrarily choose the bottom-left since they're all the same
    return
end

% Evaluate the function on the closest corners
fQ11 = f(B, L);
fQ12 = f(T, L);
fQ22 = f(T, R);
fQ21 = f(B, R);

% Coordinates of closest corners
x1 = x(L);
x2 = x(R);
y1 = y(B);
y2 = y(T);

% Check if xq is in x so that x2 - x1 = 0.
if (L == R && T ~= B)
    % Interpolate only in the y direction
    deltay = y2 - y1;
    fxy1 = fQ11; % Arbitrarily pick Q11 over Q21 since L = R
    fxy2 = fQ22; % Arbitrarily pick Q22 over Q12 since L = R
    fq = ((y2 - yq) / deltay) * fxy1 + ((yq - y1) / deltay) * fxy2;
    return
elseif (L ~= R && T == B)
    % Interpolate only in the x direction
    deltax = x2 - x1;
    fx1y = fQ11; % Arbitrarily pick Q11 over Q21 since T = B
    fx2y = fQ22; % Arbitrarily pick Q22 over Q12 since T = B
    fq = ((x2 - xq) / deltax) * fx1y + ((xq - x1) / deltax) * fx2y;
	return
end
% Default (no matches) - per Wikipedia,
fq = [x2 - xq, xq - x1]*[fQ11, fQ12; fQ21, fQ22]*[y2 - yq; yq - y1] / ((x2 - x1) * (y2 - y1));