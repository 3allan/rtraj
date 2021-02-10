function MDF = wrapMDF(func, L, resolution, spacingType, ...
    relUncertainty, randLimits, xOverLInterval)
% 
% Matt Werner (m.werner@vt.edu) - Feb 7, 2021
% 
% Provide a wrapper to clean up declarations of the mass density function
% (MDF) when defining the rocket.
% 
%    Inputs:
% 
%              func - Anonymous function that defines the nominal mass
%                     density profile.
%                     Size: 1-by-1 (function handle)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%                 L - Length of the component whose density is being
%                     described.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%        resolution - The amount of elements to include in the declaration
%                     of the sample points where the nominal mass density
%                     function will be evaluated.
%                     Size: 1-by-1 (integer)
%                     Units: - (unitless)
% 
%       spacingType - Indication as to which type of spacing to use for the
%                     sample points. Regardless of the type, the sample
%                     points (of which there are 'resolution' many), are
%                     spaced between 0 and 'L'.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                     "Linear" - Provides linear spacing between the sample
%                                points such that the distance between
%                                successive points remains constant. The
%                                sample density is, therefore, constant.
% 
%                     "Cosine" - Provides cosine spacing between the sample
%                                points such that points in the midsection
%                                near L/2 receive significant less sample
%                                points than those regions near 0 and L.
%                                This option is useful if the geometry
%                                changes rapidly near these points.
% 
%    relUncertainty - Optional(!) Relative uncertainty in the (nominal)
%                     mass density function, f.
%                     Size: 1-by-1 (scalar) OR m-by-1 (vector)
%                     Units: - (unitless)
% 
%        randLimits - Optional(!) Indication as where to treat the
%                     indicated uncertainty. Matlab random numbers are
%                     provided by default on the interval (0, 1), but this
%                     range might vary to include symmetry (-1, 1) or some
%                     sort of bias, e.g. (0.5, 0.7).
%                     Size: 1-by-2 (vector) OR m-by-2 (matrix)
%                     Units: - (unitless)
% 
%    xOverLInterval - Optional(!) Indication(s) as to where along the
%                     body the uncertainty should be applied normalized on
%                     the body's length. Therefore, entries in this
%                     quantity are interpreted as a percentage along the
%                     body's long axis.
%                     Size: m-by-2 (matrix)
%                     Units: m (meters)
% 

% Enforce that at least the first 4 inputs are given
narginchk(4, 7)

% Check if the provided spacing type is a string
checkInput(spacingType)

% No other checks - leave errors to being thrown when computing the mass
% density function in the main routines

% Determine how to space the sample points
switch upper(spacingType)
    case "LINEAR"
        x = linspace(0, L, resolution)';
    case "COSINE"
        x = cosspace(0, L, resolution)';
    otherwise
        error("Unrecognized spacing type")
end

% Evaluate the nominal density function at the sample points
f = func(x);

% Call the mass density function
if (nargin == 4)
    MDF = defineMassDensityFunction(f, x);
else
    % Ensure that all 7 inputs are given
    narginchk(7, 7)
    
    % Unnormalize the specified regions where the uncertainty should apply
    xInterval = L*xOverLInterval;
    
    MDF = defineMassDensityFunction(f, x, relUncertainty, randLimits, xInterval);
end


