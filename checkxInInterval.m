function checkxInInterval(x, I)
% 
% Matt Werner (m.werner@vt.edu) - Dec 2, 2020
% 
% Verify that the given quantity x is contained in the specified interval I.
% If the quantity x and interval I are of the same height (n), then each
% component of x is checked against the corresponding component of the
% interval I.
% 
%    Inputs:
% 
%                 x - Quantity to check whether or not all components are
%                     contained in the specified interval. Every component
%                     of x must contain a single real number.
%                     Size: n-by-1 (vector)
%                     Units: ? (corresponds to those in I)
% 
%                 I - Interval [a, b] with a <= b in which every component 
%                     of x must lie. If the height of I is 1, then every 
%                     component of x is checked against this single 
%                     interval. Otherwise, the components of x are checked 
%                     component-wise with the intervals given in I, 
%                     provided x and I are of the same height n.
%                     Size: n-by-2 (matrix)
%                     Units: ? (corresponds to those in x)
% 
%    Outputs:
% 
%                   -
% 

% Get number of elements of each x and I
numelx = numel(x);
numelI = numel(I);

% Check if a valid amount of intervals are provided against which to check x
if (2*numelx ~= numelI && numelI > 2)
    error("Value(s) cannot be compared to the given intervals.")
end

% Check if x is inside the provided intervals
for ii = 1:numelx
    xii = x(ii);
    Iii = min(ii, numelI/2);
    [a, b] = deal(I(Iii, 1), I(Iii, 2));
    if (b < a)
        error("Interval is in an invalid format.")
    elseif (xii < a || xii > b)
        if (numelx > 1)
            error("Quantity is outside of specified interval (row %1.f).", ii)
        end
        error("Quantity is outside of specified interval.")
    end
end
        