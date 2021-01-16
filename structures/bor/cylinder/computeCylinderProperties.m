function [M, CoM, IMoIatCoM, profile] = computeCylinderProperties(body)
% 
% Matt Werner (m.werner@vt.edu) - Jan 16, 2021
% 
% Obtain structural properties of a cylindrical body tube. The structural 
% properties of interest are the total mass, position of the center of
% mass, and mass moment of inertia matrix.
% 
% The center of mass is computed relative to one of the body's flat faces
% at the centerline using the (very standard) coordinate system oriented
% such that the x-axis is along the body's longitudinal axis/centerline.
% The y-axis is defined as usual and, therefore, the z-axis is oriented
% outwards.
% 
% The body's mass moment of inertia matrix is calculated relative to the
% body's center of mass (NOT the flat face at the centerline) and the
% coordinate system having the same orientation as the one used in defining
% the position of the body's center of mass, but simply shifted from the
% body's flat face at the centerline to its center of mass.
% 
% Cylinder/body profile:
%   The defining surface is given by
%                       f(x) = R,   (0 < R)
%   where R is the radius of the cylinder.
%   To make a hollow 3D shape, an inner surface 0 <= g(x) < f(x) is defined
%                     g(x) = f(x) - r,   (0 < r <= R)
%   where r is the (normal) thickness of the cylinder wall.
%  
%  y |
%    |
%  R |----------------------------------------------------  f(x)
%    |
%    |
%    |
%    |
%    |
%    |
%    |
%    |____________________________________________________|__ x  (revolve)
%  (.)                                                    L
%  z
%    Inputs:
% 
%              body - Cylinder/body object containing all relevant information
%                     used in defining the body's overall geometry and
%                     characteristics.
%                     Size: ? (structure)
%                     Units: SI
% 
% ------ Overview of content contained within the frustum ------
% 
%                OD - Outer diameter of the cylinder.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                ID - Inner diameter of the cylinder.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
% 
%            length - Length of the cylinder from face to face. For typical
%                     cylinders, this quantity is usually called the
%                     height.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%           density - Material density of the body. This density may be
%                     specified as a density function or simple scalar.
%                     Size: 1-by-1 (scalar) OR n-by-2 (matrix)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
% ---------------------------------------------------------------
% 
%   Outputs:
% 
%                 M - The cylinder's total mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               CoM - The position of the cylinder's center of mass
%                     relative to the frustum coordinate system
%                     (positioned at the front face with its x-axis pointed
%                     along the frustum's longitudinal axis towards the
%                     center of mass).
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%             IMoIG - The cylinder's mass moment of inertia matrix taken
%                     with respect to the cylinder's center of mass and the
%                     cylinder coordinate system (positioned along the
%                     longitudinal axis on the forward face with its x-axis
%                     pointed along the cylinder's longitudinal axis towards
%                     the center of mass).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%           profile - Geometric representation of the cylinder's shape
%                     according to the provided defining parameters.
%                     Size: ? (structure)
%                     Units: SI
% 

% Unpack the body components
outerRadius = body.OD/2;
innerRadius = body.ID/2;
length = body.length;
density = body.density;

%% General checks
% Check that the cylinder is physical
if (outerRadius < innerRadius || length < 0)
    error("Body (cylinder) is unphysical.")
end

% Check the density is either a scalar or an n-by-2 array with n > 1
sizedensity = size(density);
if (~(isscalar(density) || (sizedensity(1) > 1 && sizedensity(2) == 2)))
    error("Specified density (%1.0f-by-%1.0f) must be 1-by-1 or n-by-2 (where n > 1).", sizedensity(1), sizedensity(2))
end

%% Evaluate the inner/outer surfaces and density function
% Define the sampling points
numelx = 1e2; % DEFAULT value for cylinder sampling resolution
x = linspace(0, length, numelx)';

% Distribute the density appropriately
if (~isscalar(density))
    xsamp = density(:, 1);
    densitysamp = density(:, 2);
    
    % Check the provided density sampling points are well-ordered
    if (any(diff(xsamp) <= 0))
        error("Sampling points are not well-ordered.")
    end
    
    % Check that the provided density sampling points stay within the
    % bounds of the body
    checkxInInterval([xsamp(1,1), xsamp(end,1)], [0, length])
    
    % Get the number of provided elements
    numelxsamp = numel(xsamp);
    
    % Check if xsamp is more highly refined than the default sampling
    if (numelxsamp > numelx)
        % Redefine the sampling points to be more resolved than those
        % points defining the density function
        x = linspace(0, length, 2*numelxsamp);
    end
    % Override the density n-by-2 input
    density = interp1(xsamp, densitysamp, x, 'linear');
end

% Evaluate the outer and inner surfaces of the cylinder
outersurf = outerRadius*ones(size(x));
innersurf = innerRadius*ones(size(x));

%% Properties
% Compute the properties of the nosecone
[M, CoM, IMoIatCoM] = computeBORProperties(x, outersurf, innersurf, density);

% Provide the nosecone's profile (x, outersurf, innersurf, density)
profile.x = x;
profile.outersurf = outersurf;
profile.innersurf = innersurf;
profile.density = density.*ones(size(x));