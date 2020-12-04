function checkCoefficients(C, S, model)
% 
% Matt Werner (m.werner@vt.edu) - Dec 2, 2020
% 
% Analyze a set of harmonic coefficients (C, S) to determine how
% coefficients have been loaded and how many are nontrivial (different from
% zero).
% 
%    Inputs:
% 
%                 C - Cosine-series harmonic coefficients of any particular
%                     (spherical harmonic) model/correction.
%                     Size: L-by-1 (vector)*
%                     Units: - (?)
% 
%                 S - Sine-series harmonic coefficients of any particular
%                     (spherical harmonic) model/correction.
%                     Size: L-by-1 (vector)*
%                     Units: - (?)
% 
%             model - Optional(!) Indicator as to which model these
%                     coefficients (C, S) pertain to. If not included, then
%                     only a possible warning can result from this check.
%                     Including the model, however, will print lines to
%                     determine how many coefficients there are and how
%                     many are nonzero.
%                     Size: 1-by-1 (string)
%                     Units: - (N/A)
% 
%    Outputs:
% 
%                   - 
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2,
%                     where (N, M) are the maximum degree and order of the
%                     model respectively.
% 

if (~all(size(C) == size(S)))
    warning("Check model to ensure that mismatched size of coefficients is intended.")
end

if (nargin == 3)
    checkInput(model)
    L = numel(C);
    numnonzerocoeffs = numel(C(C ~= 0)) + numel(S(S ~= 0));
    disp(model + " model contains " + 2*L + " coefficients")
    disp(model + " model contains " + numnonzerocoeffs + " nonzero coefficients")
end