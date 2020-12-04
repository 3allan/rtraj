function Pnm = legendreSeries(n, m, x)
% 
% Trash
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the associated Legendre polynomial of degree n and order m
% through its closed-form series expression.
% 
%    Inputs:
% 
%                 n - Desired degree of the associated Legendre polynomial.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 
%                 m - Desired order of the associated Legendre polynomial.
%                     Size: 1-by-1 (scalar)
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

% Check that 0 <= m <= n
checkxInInterval(n, [0, inf])
if (n > 0)
    checkxInInterval(m, [0, n])
elseif (n == 0 && m ~= 0)
    error("Order m with degree n = 0 has to be zero.")
end
% Leave integer checks to those used in factorial.
% Check if x is in the proper interval
checkxInInterval(x, [-1, 1])

% % Begin recurrence relation
% if (x^2 < 0.97) % Numerically stable - upward recurrence
%     Pnn = (-1)^n * doubleFactorial(2*n - 1) * (1 - x^2
% else
%     
% end

% Try and see how fast series solution is
PnmseriesTerm = 0;
for k = m:n
    PnmseriesTerm = PnmseriesTerm + (factorial(k) / factorial(k - m)) * ...
        nchoosek(n, k) * generalizednchoosek((n + k - 1)/2, n) * ...
        x^k * x^-m;
end
Pnm = (-1)^m * 2^n * (1 - x^2)^(m/2) * PnmseriesTerm;
% ...It's fast but not very accurate