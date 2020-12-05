function checkCoefficients(C, S, model, Co, So)
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
%                Co - Optional(!) Cosine-series harmonic coefficients of
%                     any particular (spherical harmonic) model/correction
%                     corresponding with C (and S). These coefficients are
%                     the raw, uncorrected cosine coefficients of the model
%                     (but with a chosen maximum degree and order) to be 
%                     compared with C and determine how many coefficients 
%                     have been changed under the corrections.
%                     Size: L-by-1 (vector)*
%                     Units: - (?)
% 
%                So - Optional(!) Sine-series harmonic coefficients of any
%                     particular (spherical harmonic) model/correction
%                     corresponding with (C and) S. These coefficients are
%                     the raw, uncorrected sine coefficients of the model
%                     (but with a chosen maximum degree and order) to be 
%                     compared with C and determine how many coefficients 
%                     have been changed under the corrections.
%                     Size: L-by-1 (vector)*
%                     Units: - (?)
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
    L = numel(C);
    checkInput(model)
    numnonzerocoeffs = nnz(C) + nnz(S);
    disp(model + " model contains " + numnonzerocoeffs + " nonzero coefficients (" + 2*L + " total)")
elseif (nargin == 5)
    % Compare the coefficients C with Co and S with So
    % Ensure first that they are the same sizes
    checkCoefficients(C, Co)
    checkCoefficients(S, So)
    % Count the number of nonzero elements
    numnonzeroCminusCo = nnz(C - Co);
    numnonzeroSminusSo = nnz(S - So);
    numcoeffmodified = numnonzeroCminusCo + numnonzeroSminusSo;
    disp(model + " model modified " + numcoeffmodified + " coefficients")
end