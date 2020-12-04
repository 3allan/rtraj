function doublefactorialn = doubleFactorial(n)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the double factorial of the integer n
% 
%    Inputs:
%                 n - Integer that can be -1, 0, 1, 2, ...
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%    Outputs:
% 
%  doublefactorialn - Double factorial of n, namely n!!, which is NOT
%                     equivalent to (n!)! (the factorial of the factorial
%                     of n). The double factorial is defined such that for
%                     n = -1 or n = 0, one has
%                                        -1!! = 0!! = 1.
%                     Otherwise, it follows the pattern of the usual
%                     factorial but subtracts by 2 every time and stops at
%                     1 if n has odd parity or 2 if n has even parity.
%                     Thus,
%                               n!! = n (n - 2) ... 5 3 1,  (n even)
%                               n!! = n (n - 2) ... 6 4 2,  (n odd).
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 

% Check if n is an integer that is -1 or greater.
if (floor(n) ~= n || n < -1)
    error("n must be an integer greater than or equal to -1.")
end

% Handle n = -1 and n = 0
if (n == -1 || n == 0 || n == 1)
    doublefactorialn = 1;
    return
end

% Determine if n is even or odd
isodd = mod(n, 2);

% Compute the double factorial
doublefactorialn = 1;
if (isodd)
    while (n ~= 1)
        doublefactorialn = doublefactorialn * n;
        n = n - 2;
    end
else
    while (n ~= 0)
        doublefactorialn = doublefactorialn * n;
        n = n - 2;
    end
end