function MDF = defineMassDensityFunction(f, x, ...
    relUncertainty, randLimits, xInterval)
% 
% Matt Werner (m.werner@vt.edu) - Feb 6, 2021
% 
% Define the (M)ass (D)ensity (F)unction for a component of the vehicle's
% body. This quantity, MDF, is to be used in determining the component's
% structural properties (mass, center of mass, and mass moment of inertia
% matrix).
% 
% Permitting the density to be specified as a function allows for
% possibilities of including dynamical effects of inhomogeneous materials.
% 
% Examples of applications: a nosecone's tip may be made from a different
% material from the rest of its body, a nosecone may be carrying a payload,
% a casing may be carrying internal components (electronics, nozzle, etc.),
% uncertainties may exist in the body's material/payload/other internal
% components.
% 
%    Inputs:
% 
%                 f - Nominal mass density function evaluated at the sample
%                     points 'x'. Further adjustments may be made to this
%                     quantity through the use of uncertainties to
%                     accomodate for unknown/internal masses.
%                     Size: n-by-1 (vector)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%                 x - Sample points at which to evaluate the nominal mass
%                     density function.
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
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
%         xInterval - Optional(!) Indication(s) as to where along the
%                     body the uncertainty should be applied.
%                     Size: m-by-2 (matrix)
%                     Units: m (meters)
% 
%    Outputs:
% 
%               MDF - The body's mass density function which includes
%                     any effects due to payloads and uncertainty.
%                     Size: n-by-1 (vector)
%                     Units: kg/m3 (kilograms per cubic meter)
% 

% Check first if there is only 1 input and it makes sense
if (nargin >= 1 && isscalar(f) && f > 0)
    MDF = f;
    if (nargin == 1)
        % Default behavior is simply to return the provided constant
        % density
        return
    elseif (nargin == 2)
        % Default behavior is to return the provided constant density in a
        % size that matches that of x
        MDF = MDF * ones(size(x));
        MDF = [x, MDF];
        return
    end
end

% Check that x is well-ordered
if (any(diff(x) <= 0))
    error("Invalid positions.")
end
% Check that x, f have the same dimensions
if (~all(size(x) == size(f)) && ~isscalar(f))
    error("Incompatible dimensions for mass density function.")
end


% Assign the nominal MDF
density = f;

% Perturb this result if the relative uncertainty is provided for specific
% regions within the body
if (nargin >= 3)
    % Enforce that all inputs regarding uncertainty are given
    narginchk(5, 5)
    
    % Check if the relative uncertainty is valid
    if (any(relUncertainty < 0))
        error("Invalid uncertainty.")
    elseif (all(relUncertainty == 0))
        % Perfect values
        return
    end
    
    % Allocate memory for the relative uncertainty result from this
    % procedure
    ddensity_rel = zeros(size(x));
    
    % Get the number of regions requesting an uncertainty
    numRegions = size(xInterval, 1);
    
    % Determine if the indicated relative uncertainties and the limits for
    % the random numbers only have 1 row. In such a case, assume that the
    % values are to be held constant for each region
    if (isscalar(relUncertainty))
        relUncertainty = relUncertainty * ones(numRegions, 1);
    end
    if (size(randLimits, 1) == 1)
        randLimits = randLimits .* ones(numRegions, 2);
    end
    
    % Begin adding uncertainties to the indicated regions
    for region = 1:numRegions    
        % Check that the interval is valid
        random_lower = randLimits(region, 1);
        random_upper = randLimits(region, 2);
        if (random_lower > random_upper || abs(random_lower) > 1 || abs(random_upper) > 1)
            error("Invalid interval.")
        end
        % Update the default interval measures and offset, (a, b)
        intervalMeasure = random_upper - random_lower;
        intervalOffset = random_lower;
    
        % Determine where the uncertainties should be applied along the
        % body
        % 
        % Check that the interval is valid
        checkxInInterval(xInterval(region, 1:2)', [0, max(x)]);
        if (diff(xInterval(region, 1:2)) < 0)
            error("Invalid error (row %1.0f).", region)
        end
        % Unpack the interval
        xUncertainty_lower = xInterval(region, 1);
        xUncertainty_upper = xInterval(region, 2);
        % Find the limits as to where to apply the uncertainty in this
        % region
        [~, indUncertainty_lower] = min(abs(x - xUncertainty_lower));
        [~, indUncertainty_upper] = min(abs(x - xUncertainty_upper));
        region_limits = indUncertainty_lower:indUncertainty_upper;
    
        % Create a random set of absolute uncertainties along the provided
        % interval
        randomSample = rand(size(x(region_limits)));
        randomNumsOnAdjustedInterval = intervalMeasure*randomSample + intervalOffset;
        ddensity_rel(region_limits) = relUncertainty(region)*randomNumsOnAdjustedInterval;
    end
    
    % Calculate the perturbed mass density function
    MDF = [x, density.*(1 + ddensity_rel)];
end