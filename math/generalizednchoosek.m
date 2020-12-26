function xchoosek = generalizednchoosek(x, k)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the generalized binomial coefficient such that x is allowed to be
% any complex number according to
%   /   \
%   | x |      x (x - 1) (x - 2) ... (x - (k - 1))
%   |   |  =  -------------------------------------,
%   | k |                      k!
%   \   /
% where k is still a natural number (0, 1, 2, ...).
% 
%    Inputs:
% 
%                 x - Arbitrary complex number to insert into the
%                     generalized binomial formula.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 
%                 k - Natural number.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 
%    Outputs:
% 
%          xchoosek - Generalized binomial coefficient calculated using the
%                     arbitrary complex number x and natural number k.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 

% Check if k is a natural number
if (k < 0 || floor(k) ~= k)
    error("The second input has to be a non-negative integer.")
end

% Calculate the generalized binomial coefficient
xchoosek = 1;
for p =1:k
    xchoosek = xchoosek * (x + 1 - p);
end
xchoosek = xchoosek / factorial(k);