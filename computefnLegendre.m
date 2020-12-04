function Pbarnm = computefnLegendre(n, x, m)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the fully normalized associated Legendre polynomials according to
% the  standard normalization convention,
%                _______________________________          
%   _           /  (n - m)!                               
%   P  (x) =   /  ----------  (2 - d  ) (2n + 1)  P  (x), 
%    nm      \/    (n + m)!         0m             nm
% 
% provided that 0 <= m <= n and |x| <= 1 is real and where d0m is the
% Kronecker delta function, kroneckerDelta(0, m). There are forms of the
% Legendre polynomials with both subscripted n and superscripted m. The
% mixed flavor is not used in any computations here, but for clarity, 
% their relation is
%                m    m
%   P  (x) = (-1)   P  (x).
%    nm              n
% The (mixed) associated Legendre polynomial comes from the Legendre
% polynomial (of order 0)
%                                    m
%     m          m     2    m/2     d 
%   P  (x) = (-1)   (-x + 1)      -------  P (x),
%    n                            dx...dx   n
% which finally satisfy Legendre's (differential) equation.
%        2   
%  (1 - x ) P''(x) - 2x P'(x) + n(n + 1) P (x) = 0
%            n           n                n
% The procedure in solving Legendre's equation may be skipped completely by
% introducing the Rodrigues' formula, which gives the solution as a
% procedural definition involving derivatives of a 2n-th degree polynomial.
%                          n          
%                 1       d       2     n
%      P  (x) = ------  ------- (x  - 1) 
%       n       (2n)!!  dx...dx      
% The symbol (2n)!! here is the double factorial, NOT ((2n)!)!. In fact,
% one has that the double factorial is rewritten (2n)!! = 2^n * n!. So
% under this framework, the associated Legendre polynomial of degree n and
% order m, where 0 <= m <= n and -1 <= x <= 1, is 
%                                     n + m          
%                 1        2    m/2  d        2     n
%      P  (x) = ------  (-x + 1)    ------- (x  - 1) .
%       nm      (2n)!!              dx...dx       
% Note that multiple definitions of the associated Legendre polynomial Pnm
% exist, most notably the one that includes the Condon-Shortley phase
% (-1)^m, which does affect the Legendre's equation (which is why so much
% excessive detail was included here). Also note that the Legendre
% polynomial Pnm(x) is not actually a polynomial if m has odd parity.
% 
%    Inputs:
% 
%                 n - Degree of the associated Legendre polynomial.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%                 x - Value at which to evaluate the Legendre polynomials.
%                     Must satisfy the requirements that x is real and in
%                     the interval (closed) interval [-1, 1].
%                     Size: 1-by-N (vector)
%                     Units: - (unitless)
% 
%                 m - Optional(!) Order of the associated Legendre
%                     polynomial. If left unspecified, then all orders from
%                     0 to n are returned as a column vector for each input
%                     x.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 

% Check input
checkxInInterval(x, [-1, 1])
% Check if the degree is negative or not an integer
if (n < 0 || floor(n) ~= n)
    error("Order n must be positive scalar integer.")
end
% Check if order exceeds the degree
if (nargin == 3 && m > n)
    error("Order m exceeds degree n.")
end

% Compute fully normalized associated Legendre polynomial of degree n and
% order m.
Pbarnms = sqrt(2*n + 1) * legendre(n, x, 'sch');

if (nargin == 2)
    Pbarnm = Pbarnms;
    return
end
Pbarnm = Pbarnms(1 + m, :);