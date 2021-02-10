function mixedf = mixf(x, f)
% 
% Matt Werner (m.werner@vt.edu) - Dec 9, 2020
% 
% Calculate the mixed value of f according to the rule of mixtures
% 
% mixedf = x*f  + (1 - x)*f ,
%             2            1
% where x is some (real) value between 0 and 1 (inclusive). Doing so
% provides a smooth, linear way to vary between f1 and f2 depending on the
% values x.
% 
%    Inputs:
% 
%                 x - Provides a means to determine how much mixing occurs.
%                     This quantity must satisfy the requirement that it be
%                     between 0 and 1 (inclusive) for the rule of mixtures
%                     to make sense. If x is zero, then no mixing occurs,
%                     so the mixed quantity is simply the unmixed quantity
%                     (f1). If x is unity otherwise, then the mixed
%                     quantity is fully mixed and correspondingly returns
%                     the fully mixed value (f2).
%                     Size: n-by-1 (vector)
%                     Units: - (N/A)
% 
%                 f - Quantity to mix between the limiting values it
%                     contains (f1 and f2) by the amount x. The first value
%                     contained in f is the unmixed (original) quantity,
%                     and the second value contained in f is the fully
%                     mixed (changed) quantity. These two values only need
%                     be real - they may be identical, though in this case,
%                     no mixing occurs.
%                     Size: n-by-2 (matrix)
%                     Units: ?
% 

% Check if x is in the unit interval
checkxInInterval(x, [0, 1])

% Check if no mixing occurs
if (all(f(:, 1) == f(:, 2)))
    mixedf = f;
    return
end

% Calculate the mixed f
f1 = f(:, 1);
f2 = f(:, 2);
mixedf = x.*f2 + (1 - x).*f1;

