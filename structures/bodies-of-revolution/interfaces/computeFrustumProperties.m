function [profile, M, CoM, IMoIatCoM] = ...
    computeFrustumProperties(frustum, stage)
% 
% Matt Werner (m.werner@vt.edu) - Jan 16, 2021
% 
% Obtain structural properties of a frustum. The structural properties of
% interest are the total mass, position of the center of mass, and mass
% moment of inertia matrix.
% 
% The center of mass is computed relative to the frustum's forward face at
% the centerline using the (very standard) coordinate system oriented such
% that the x-axis is along the frustum's longitudinal axis/centerline. The
% y-axis is defined as usual and, therefore, the z-axis is oriented
% outwards.
% 
% The frustum's mass moment of inertia matrix is calculated relative to
% the frustum's center of mass (NOT the forward face at the centerline) and
% the coordinate system having the same orientation as the one used in
% defining the position of the frustum's center of mass, but simply shifted
% from the frustums's forward face at the centerline to its center of mass.
% 
% Frustum profile:
%   The defining surface is given by
%               f(x) = r + (R - r) * (x/L),   (0 < r < R)
%   where r is the front radius and R is the rear radius. 
%   To make a hollow 3D shape, an inner surface 0 <= g(x) < f(x) is defined
%   according to the thickness type ('thicknessType') and thickness
%   ('thickness').
% 
%   y |
%     | 
%   R |                            ,,/  f(x)
%     |                         ,,/
%     |                      ,,/
%     |                   ,,/
%     |                ,,/
%     |             ,,/
%     |          ,,/
%     |       ,,/
%     |    ,,/
%     | ,,/
%   r |/
%     |
%     |__________________|___________|________________________ x (revolve)
% z (.)                 CoM          L                          
% 
%    Inputs:
% 
%           frustum - Frustum object containing all relevant information
%                     used in defining the frustum's overall geometry and
%                     characteristics.
%                     Size: ? (structure)
%                     Units: SI
% 
% ------ Overview of content contained within the frustum ------
% 
%           frontOD - Diameter of the frustum's front (base/smaller) face.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%            rearOD - Diameter of the frustum's rear (base/larger) face.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
% 
%            length - Length of the frustum from front to rear faces. For
%                     typical frustums, this quantity is usually called the
%                     height.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%     thicknessType - Specifies how to interpret the thickness.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                   Possible options are:
%                 "Horizontal" - Thickness specifies the horizontal
%                                distance between the inside and outside
%                                surfaces. The inside surface is exactly
%                                the same as the outside surface but
%                                shifted over by 'thickness' units.
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
%                     "Normal" - Thickness specifies the distance between
%                                the outer and inner surfaces as measured
%                                relative to the direction normal to the
%                                outer surface. In this case, the thickness
%                                may only be scalar valued.
% 
%         thickness - Distance separating the frustum's outer surface from
%                     its inner surface as specified by the thickness type.
%                     Size: 1-by-1 (scalar) OR n-by-2 (vector, "pointwise")
%                     Units: m (meters)
% 
%           density - Material density of the frustum. This density may be
%                     specified as a density function or simple scalar, so
%                     long as it acts in accordance with the permissions of
%                     the thicknessType.
%                     Size: 1-by-1 (scalar) OR n-by-2 (matrix)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
% ---------------------------------------------------------------
% 
%             stage - Indication as to which stage the frustum belongs. The
%                     first stage (stage = 1) indicates the frustum that
%                     mounts to the booster that initially fires and brings
%                     the vehicle off of the launch rail. The final stage
%                     indicates the frustum that mounts the nosecone to the
%                     last sustainer stage (if such a frustum is
%                     specified).
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%   Outputs:
% 
%                 M - The frustum's total mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               CoM - The position of the frustum's center of mass
%                     relative to the frustum coordinate system
%                     (positioned at the front face with its x-axis pointed
%                     along the frustum's longitudinal axis towards the
%                     center of mass).
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%             IMoIG - The frustum's mass moment of inertia matrix taken
%                     with respect to the frustum's center of mass and the
%                     frustum coordinate system (positioned along the
%                     longitudinal axis on the forward face with its x-axis
%                     pointed along the frustum's longitudinal axis towards
%                     the center of mass).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%           profile - Geometric representation of the frustum's shape
%                     according to the provided defining parameters and
%                     requests.
%                     Size: ? (structure)
%                     Units: SI
% 

% Unpack the frustum for this stage
frontRadius = frustum.frontOD(stage, 1)/2;
rearRadius = frustum.rearOD(stage, 1)/2;
length = frustum.length(stage, 1);
thicknessType = frustum.thicknessType(stage, 1);
thickness = frustum.thickness(stage, 1);
density = frustum.density(stage, 1);

% Check if all values are 0 or NaN
tmp_els = [frontRadius, rearRadius, length, thickness, density];
if (all(tmp_els == 0) || all(isnan(tmp_els)))
    % No frustum - return 0 for the outputs
    [M, CoM, IMoIatCoM, profile] = deal(0);
    return
end

% Check the thicknessType
checkInput(thicknessType)

% Check the thickness
checkxInInterval(thickness, [0, inf])

% Check that the rear and front faces have been filled out correctly
if (frontRadius > rearRadius)
    error("Backwards frustum.")
end
if (frontRadius == rearRadius)
    error("Frustum may not be a cylinder.")
end

% Check the density is either a scalar or an n-by-2 array with n > 1
sizedensity = size(density);
if (~(isscalar(density) || (sizedensity(1) > 1 && sizedensity(2) == 2)))
    error("Specified density (%1.0f-by-%1.0f) must be 1-by-1 or n-by-2 (where n > 1).", sizedensity(1), sizedensity(2))
end

%% Outer Surface
% Define the sampling points
x = linspace(0, length, 1e2)';

% Define the frustum profile
fun = @(x) frontRadius + (rearRadius - frontRadius)./length.*x;

% Evaluate the outer surface of the frustum
outersurf = fun(x);

%% Inner Surface
% Determine the inner surface appropriately according to the
% requested thickness type
switch upper(thicknessType)
    case "HORIZONTAL"
        % Check that the thickness will result in a frustum with two open
        % ends
        if (thickness > rearRadius*length / (rearRadius + frontRadius))
            error("Frustum is closed at its forward face")
        end
        % Shift the outside surface by 'thickness' horizontally
        innersurf = fun(x - thickness);
        
    case "VERTICAL"
        % Check that the thickness will result in a frustum with two open
        % ends
        if (thickness >= frontRadius)
            error("Frustum is closed at its forward face")
        end
        % Shift the outside surface by 'thickness' vertically
        innersurf = outersurf - thickness;
        
    case "VERTICAL-POINTWISE"
        % Apply a vertical shift downwards from the outer surface
        % by an amount specified by the n-by-2 array 'thickness',
        % where the first column corresponds with the physical
        % location within the frustum at which to apply the
        % thickness indicated in the second column. The sizes must
        % match, so update the outer surface to be more fine and
        % (linearly) interpolate to match sizes.
        maxElements = max(numel(x), size(thickness, 1));
        xintp = linspace(0, length, 100 + maxElements)';
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
        % Find the angle of the frustum
        heightGain = rearRadius - frontRadius;
        sideLength = sqrt(heightGain^2 + length^2);
        sintc = heightGain/sideLength;
        
        % Calculate the vertical thickness
        verticalThickness = thickness/sintc;
        
        % Shift the outside surface by the equivalent vertical thickness
        innersurf = outersurf - verticalThickness;
        
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