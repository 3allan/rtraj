function PnmRecurrence = computeLegendreRecurrence(n, m, x)
% 
% Trash
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the recurrence relations for the associated Legendre polynomials
% of degree n and order m (given) to advance to degree n + 1 having maximum
% order m + 1.
% 
%    Inputs:
% 
%               Pnm - Legendre polynomials of degree n and orders 0 to n
%                     evaluated at the specified x.
%                     Size: (n+1)-by-1 (vector)
%                     Units: - (N/A)
% 
%                 x - Argument at which the associated Legendre polynomials
%                     Pnm are evaluated.
%                     Size: 1-by-1 (scalar)
%                     Units: ?
% 
%    Outputs:
% 
%     PnmRecurrence - Legendre polynomials of degree n + 1 and orders 0 to
%                     n + 1 evaluated at the specified x.
%                     Size: (n+2)-by-1 (vector)
%                     Units: - (N/A)
% 

% Check that 0 < m < n
checkxInInterval(n, [0, inf])
checkxInInterval(m, [0, n])
% Check if x is in the proper interval
checkxInInterval(x, [-1, 1])

% Try to do a recurrance
if (n == 0 && m == 0)
    PnmRecurrence = 1;
elseif (n == 1 && m == 0)
    PnmRecurrence = x;
elseif (n == 2 && m == 0)
    PnmRecurrence = 1.5*x^2 - 0.5;
elseif (m == n)
    PnmRecurrence = (-1)^n * doubleFactorial(2*n - 1) * (1 - x^2)^(0.5*n);
elseif (m == 0)
    PnmRecurrence = ((2*n - 3)*x*computeLegendreRecurrence(n-2, 0, x) - (n - 2)*computeLegendreRecurrence(n-3, 0, x)) / (n - 1);
else
    PnmRecurrence = -0.5*sqrt(1-x^2)/m * (computeLegendreRecurrence(n-1, m+1, x) + (n + m - 1)*(n + m)*computeLegendreRecurrence(n-1, m-1, x));
end

