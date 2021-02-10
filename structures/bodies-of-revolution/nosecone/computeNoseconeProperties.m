function [profile, M, CoM, IMoIatCoM] = computeNoseconeProperties(nosecone)
% 
% Matt Werner (m.werner@vt.edu) - Jan 15, 2021
% 
% Obtain structural properties of a nosecone. The structural properties of
% interest are the total mass, position of the center of mass, and mass
% moment of inertia matrix.
% 
% The center of mass is computed relative to the nosecone's tip using the
% (very standard) coordinate system oriented such that the x-axis is along
% the nosecone's longitudinal axis/centerline. The y-axis is defined as
% usual and, therefore, the z-axis is oriented outwards.
% 
% The nosecone's mass moment of inertia matrix is calculated relative to
% the nosecone's center of mass (NOT the vertex) and the coordinate system
% having the same orientation as the one used in defining the position of
% the nosecone's center of mass, but simply shifted from the nosecone's
% vertex to its center of mass.
% 
% Nosecone profile:
%   The defining surfaces f(x) are given according to the nosecone series.
%   To make a hollow 3D shape, an inner surface 0 <= g(x) < f(x) is defined
%   according to the thickness type ('thicknessType') and thickness
%   ('thickness').
% 
%   y |
%     |                                            
%   R |                                     ,_____________  f(x)
%     |                              ,_____/
%     |                        ,____/
%     |                   ,___/
%     |               ,__/
%     |           ,__/
%     |        ,_/
%     |     ,_/
%     |   ,/
%     | ,/               
%   r |/______________________________|__________________|_____ x (revolve)
% z (.)                              CoM                 L     
% 
%    Inputs:
% 
%          nosecone - Nosecone object containing all relevant information
%                     used in defining the nosecone's overall geometry and
%                     characteristics.
%                     Size: ? (structure)
%                     Units: SI
% 
% ------ Overview of content contained within the nosecone ------
% 
%            radius - Radius of the nosecone's base (a circle) as measured
%                     from the nosecone's longitudinal axis/centerline.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%            length - Length of the cone from its tip to its base. For
%                     typical (right, circular) cones, this quantity is
%                     usually called the height.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%     thicknessType - Specifies how to interpret the thickness.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                   Possible options are:
%                    "Frontal" - Thickness specifies the distance between
%                                the inside and outside vertices. The
%                                inside surface is exactly the same as the
%                                outside surface but shifted over by
%                                'thickness' units. Nosecones are typically
%                                convex shapes, so unphysical shapes are
%                                usually difficult to generate (unphysical
%                                shapes occur if the inside surface
%                                intersects the outside surface,
%                                discontinuous, etc.)
% 
%                   "Vertical" - Thickness specifies the "vertical" (in the
%                                y-direction as defined by the BOR
%                                coordinate system) separation between the
%                                inner and outer surfaces. The inner
%                                surface is simply the outer surface minus
%                                the specified thickness (scalar). Any
%                                inner surface falling below the cone's
%                                longitudinal axis is ignored such that the
%                                region between the outer vertex and inner
%                                vertex is solid.
% 
%         "Vertical-Pointwise" - Same description as for "Vertical" except
%                                the thickness is specified as a n-by-2
%                                array, where the first column consists of
%                                the sampling points at which the thickness
%                                applies (x) and the second column is the
%                                thickness. The thickness is linearly 
%                                interpolated between sampling points at
%                                the sampling points of the outer surface.
% 
%                     "Normal" - Not yet implemented.
% 
%         thickness - Distance separating the nosecone's outer surface from
%                     its inner surface as specified by the thickness type.
%                     Size: 1-by-1 (scalar) OR n-by-2 (vector, "pointwise")
%                     Units: m (meters)
% 
%           density - Material density of the nosecone. A scalar density
%                     indicates that the density applies everywhere (at all
%                     points) along the nosecone. Otherwise, the density
%                     may be specified as a n-by-2 array indicating the
%                     material density (column 2) at specified sampling
%                     points along the nosecone (column 1). In this case,
%                     the density is linearly interpolated to be evaluated
%                     at the sampling points defining the outer surface.
%                     Size: 1-by-1 (scalar) OR n-by-2 (matrix)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%            series - Cone classification.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                     "Linear" - Indicates a right circular cone having a
%                                linear (completely straight) profile with
%                                an infinitely sharp vertex. The linear
%                                profile is precisely the same as a
%                                parabolic profile with its cone constant
%                                set to 0.
% 
%                      "Power" - Indicates a cone having a power profile
%                                such that its profile is curved as found
%                                in functions of positive fractional powers
%                                less than unity (1) and is parameterizable
%                                by a cone constant k.
% 
%                  "Parabolic" - Indicates a cone having a parabolic
%                                (quadratic) profile such that its profile
%                                is curved as a parabola and is
%                                parameterizable by a cone constant k.
% 
%              "Tangent-Ogive" - Indicates a cone having a circular profile
%                                such that its profile is taken as the arc
%                                from a circle defined by the nosecone's
%                                length and radius. This nosecone is
%                                unparameterizable.
% 
%               "Secant-Ogive" - Indicates a cone having a circular profile
%                                such that its profile is taken as the arc
%                                from a circle defined to be greater in
%                                radius than the circle defining the
%                                tangent-ogive for nosecones of the same
%                                length and radius. This nosecone is
%                                parameterizable.
% 
%                      "Haack" - Indicates a cone having a mathematically
%                                derived profile for minimizing aerodynamic
%                                drag and is parameterizable by a cone
%                                constant k. This series includes the von
%                                Karman (vK) ogive.
% 
%                 k - Cone constant determining the cone's curvature.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%                    Extended description for specific series:
%                     "Linear" - N/A
% 
%                      "Power" - The cone constant is restricted to be
%                                within the (half-closed) interval
%                                              0 < k <= 1.
%                                Having k = 0 results in a cylinder, which
%                                is why 0 is not included in the interval.
%                                At k = 1/2, the nosecone is a paraboloid.
%                                At k = 1, the nosecone reduces back to the
%                                right circular cone of the linear series
%                                (which is also a special case of the
%                                parabolic series).
% 
%                  "Parabolic" - The cone constant k is restricted to be
%                                within the (closed) interval
%                                              0 <= k <= 1.
%                                Having k = 0 results in a shape
%                                identifying precisely with a right
%                                circular cone having an infinitely sharp
%                                vertex. Increasing the cone constant
%                                upwards towards its upper limit of 1
%                                increases the cone's curvature. In this
%                                case, the vertex is still sharp, but the
%                                sides at least curve smoothly into it. 
%                                At k = 1, the cone fully curves over such
%                                that the profile is parallel with its
%                                longitudinal axis at the cone's base.
% 
%              "Tangent-Ogive" - N/A
% 
%               "Secant-Ogive" - The cone constant k is restricted to the
%                                open interval (, inf)
%                                     (R2 + L2)/2R < k < infinity,
%                                where R is the nosecone's base radius and
%                                L is the nosecone's overall length from
%                                tip to base. This open lower limit is
%                                enforced since, upon equality, the
%                                nosecone would actually be a tangent
%                                ogive.
% 
%                      "Haack" - The cone constant k is technically
%                                restricted within the (half-closed)
%                                interval [0, inf)
%                                           0 <= k < infinity,
%                                though cone constants (typically) greater
%                                than 1 give impractical nosecones.
%                                Having k = 0 results in the von Karman
%                                (vK) ogive cone. Increasing the cone
%                                constant increases the cone's curvature.
%                                At k = 2/3, the Haack series nosecone
%                                fully curves over such that its profile is
%                                parallel with its longitudinal axis at the
%                                cone's base.
%                                For k > 2/3, the profile will continue to
%                                bulge outwards, exceeding the base radius
%                                between the tip and base.
% 
% ---------------------------------------------------------------
% 
%   Outputs:
% 
%                 M - The nosecone's total mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               CoM - The position of the nosecone's center of mass
%                     relative to the nosecone coordinate system
%                     (positioned at the tip with its x-axis pointed along
%                     the nosecone's longitudinal axis towards the center
%                     of mass).
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%         IMoIatCoM - The nosecone's mass moment of inertia matrix taken
%                     with respect to the nosecone's center of mass and the
%                     nosecone coordinate system (positioned at the tip
%                     with its x-axis pointed along the conecone's
%                     longitudinal axis towards the center of mass).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%           profile - Geometric representation of the nosecone's shape
%                     according to the provided defining parameters and
%                     requests.
%                     Size: ? (structure)
%                     Units: SI
% 

% Unpack elements from the nosecone structure
radius = nosecone.OD/2;
length = nosecone.length;
thicknessType = nosecone.thicknessType;
thickness = nosecone.thickness;
density = nosecone.density;
series = nosecone.series;
k = nosecone.k;

% Check the series and thicknessType to ensure that they are strings
checkInput(series)
checkInput(thicknessType)

% Check the radius, length, and thickness for positive-definiteness
checkxInInterval([radius; length; thickness], [0, inf])

% Check the density is either a scalar or an n-by-2 array with n > 1
sizedensity = size(density);
if (~(isscalar(density) || (sizedensity(1) > 1 && sizedensity(2) == 2)))
    error("Specified density (%1.0f-by-%1.0f) must be 1-by-1 or n-by-2 (where n > 1).", sizedensity(1), sizedensity(2))
end

%% Outer surface

% Define the sampling points by first defining angular spacing coordinates
x = cosspace(0, length, 1e3)';

% Define the structural properties of the nosecone based on its type
% (series) and how the thickness is applied
switch upper(series)
    case "LINEAR"
        % Define the linear profile
        fun = @(x) (radius/length)*x;
    case "POWER"
        % Check the power cone constant
        checkxInInterval(k, [0, 1])
        if (k == 0), error("Nosecone geometry is invalid."), end
        
        % Define the power profile
        fun = @(x) radius * (x/length).^k;
    case "PARABOLIC"
        % Check the parabolic cone constant
        checkxInInterval(k, [0, 1])
        
        % Define the parabolic profile
        fun = @(x) radius * (2 - k * x/length) .* x/length / (2 - k);
    case "TANGENT-OGIVE"
        % Calculate the defining circle's radius
        rCirc = (radius^2 + length^2)/(2*radius);
        
        % Define the tangent ogive profile
        fun = @(x) radius + sqrt(rCirc^2 - (length - x).^2) - rCirc;
    case "SECANT-OGIVE"
        % Recase the cone constant to be the radius of the secant ogive's
        % defining circle
        rCirc = k;
        
        % Check that the cone constant (the radius of the defining circle)
        % is valid. The cone constant is valid if the resulting defining
        % circle exhibits a radius greater than the defining circle for a
        % tangent ogive of equal nosecone radii and overall lengths
        radius2plusLength2 = radius^2 + length^2;
        if (rCirc <= radius2plusLength2/(2*radius))
            error("Invalid cone constant for secant ogive.")
        end
        
        % Calculate the auxillary variable
        a = acos(sqrt(radius2plusLength2)/(2*rCirc)) - atan(radius/length);
        
        % Define the secant ogive profile
        fun = @(x) sqrt(rCirc^2 - (rCirc*cos(a) - x).^2) - rCirc*sin(a);
    case "HAACK"
        % Provide an auxillary anonymous function that undoes the
        % transformation from t to x. The Haack series will be the only
        % series that references this anonymous function from the main
        % profile function 'fun'. Undoing this transformation is convenient
        % so that every series may continue using 'fun' without having to
        % provide a complicated formula for the Haack series
        Haacktheta = @(x) acos(1 - 2*x/length);
        
        % Define the Haack series profile
        fun = @(x) radius * sqrt((1/pi) * (Haacktheta(x) - ...
            sin(2*Haacktheta(x))/2 + k*sin(Haacktheta(x)).^3));
    otherwise
        error("Nosecone series not recognized.")
end

% Evaluate the outer surface of the cone with this profile
outersurf = fun(x);
        
%% Inner surface
% Determine the inner surface appropriately according to the
% requested thickness type
switch upper(thicknessType)
    case "FRONTAL"
        % Shift the outside surface by 'thickness' and ignore any
        % points falling beneath the x-axis (become negative
        % height) by applying the Heaviside function. At x = 0, the
        % profile has fun(0) = 0, so Matlab's interpretation of the
        % Heaviside function having H(0) = 0.5 is irrelevant.
        xminust = x - thickness;
        innersurf = fun(xminust).*heaviside(xminust);
    case "VERTICAL"
        % Apply a vertical shift downwards from the outer surface
        % by an amount 'thickness'
        innersurf = outersurf - thickness;
        innersurf(innersurf < 0) =  0;
    case "VERTICAL-POINTWISE"
        % Apply a vertical shift downwards from the outer surface
        % by an amount specified by the n-by-2 array 'thickness',
        % where the first column corresponds with the physical
        % location within the nosecone at which to apply the
        % thickness indicated in the second column. The sizes must
        % match, so update the outer surface to be more fine and
        % (linearly) interpolate to match sizes.
        maxElements = max(numel(x), size(thickness, 1));
        xintp = linspace(0, length, 500 + maxElements)';
        outersurf = fun(xintp);

        % Unpack the thickness array
        xsamp = thickness(:, 1);
        thicknesssamp = thickness(:, 2);

        % Interpolate to match the new size of x
        thicknessintp = interp1(xsamp, thicknesssamp, xintp, 'linear');

        % Apply the pointwise vertical thickness specification
        innersurf = outersurf - thicknessintp;
        
        % Update the position vector x to be xintp to keep all sizes the
        % same
        x = xintp;
    case "NORMAL"
        error("Normal thickness not yet implemented.")
    otherwise
        error("Thickness type not recognized.")
end

%% Density
% Distribute the density appropriately
if (~isscalar(density))
    xsamp = density(:, 1);
    densitysamp = density(:, 2);
    % Override the density n-by-2 input
    density = interp1(xsamp, densitysamp, x, 'linear');
end

%% Properties
% Compute the properties of the nosecone
[M, CoM, IMoIatCoM] = computeBORProperties(x, outersurf, innersurf, density);

% Provide the nosecone's profile (x, outersurf, innersurf, density)
profile.x = x;
profile.outersurf = outersurf;
profile.innersurf = innersurf;
profile.density = density.*ones(size(x));