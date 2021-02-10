function [profile, M, CoM, IMoIatCoM, dIMoIatCoMdt] = ...
    computeMotorProperties(motor, stage)
% 
% Matt Werner (m.werner@vt.edu) - Feb 8, 2021
% 
% Obtain structural properties of a cylindrical motor. The structural
% properties of interest are the mass, center of mass, and mass moment of
% inertia. The motor's grain geometry is assumed to be of the BATES
% variety, which is a hollow circular cylinder.
% 
% The center of mass is computed relative to one of the body's flat faces
% at the centerline using the (very standard) coordinate system oriented
% such that the x-axis is along the motor's longitudinal axis/centerline.
% The y-axis is defined as usual and, therefore, the z-axis is oriented
% outwards.
% 
% The motor's mass moment of inertia matrix is calculated relative to the
% motor's center of mass (NOT the flat face at the centerline) and the
% coordinate system having the same orientation as the one used in defining
% the position of the motor's center of mass, but simply shifted from the
% motor's flat face at the centerline to its center of mass.
% 
% Motor profile:
%   The defining surface is given by
%                       f(x) = R,   (R > 0)
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
%             motor - Motor object containing all relevant information
%                     used in defining the motor's overall geometry and
%                     characteristics.
%                     Size: ? (structure)
%                     Units: SI
% 
% -------- Overview of content contained within the motor --------
% 
%                OD - Outer-diameter of the motor.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%             IDatt - Current inner-diameter of the motor. The motor is
%                     assumed to be a BATES grain, so the inner diameter is
%                     constant in space but not in time.
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
%           massatt - Current mass of the motor.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
% ---------------------------------------------------------------
% 
%             stage - Indication as to which stage the motor belongs. 
%                     The first stage (stage = 1) indicates the motor
%                     belonging to the booster that initially fires and
%                     brings the vehicle off of the launch rail. The final
%                     stage indicates the motor within the cylinder/casing
%                     that mounts to the nosecone and is the last sustainer
%                     to fire.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%   Outputs:
% 
%                 M - The motor's total mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               CoM - The position of the motor's center of mass relative
%                     to the motor coordinate system (positioned at the
%                     front face with its x-axis pointed along the motor's
%                     longitudinal axis towards the center of mass).
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%         IMoIatCoM - The motor's mass moment of inertia matrix taken
%                     with respect to the motor's center of mass and the
%                     motor coordinate system (positioned along the
%                     longitudinal axis on the forward face with its x-axis
%                     pointed along the motor's longitudinal axis towards
%                     the center of mass).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%      dIMoIatCoMdt - The motor's mass moment of inertia matrix taken
%                     with respect to the motor's center of mass and the
%                     motor coordinate system (positioned along the
%                     longitudinal axis on the forward face with its x-axis
%                     pointed along the motor's longitudinal axis towards
%                     the center of mass).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%           profile - Geometric representation of the motor's shape
%                     according to the provided defining parameters.
%                     Size: ? (structure)
%                     Units: SI
% 

% Unpack the motor components
outerRadius = motor.OD(stage, 1)/2;
innerRadius = motor.ID(stage, 1)/2;
innerRadiusatt = motor.IDatt(stage, 1)/2;
length = motor.length(stage, 1);
mass = motor.mass(stage, 1);
burnTime = motor.burnTime(stage, 1);

% The motor mass is specified rather than its mass density
M = mass;

% % Assume that the density of the grain is constant such that mass/volume
% is constant. In this case, the center of mass stays constant at
% midlength.
CoM = [length/2; 0; 0];

% With a constant density assumption, then the inertia matrix about the
% coordinate system's origin shown in the description is the following
sumOfSquares = innerRadiusatt^2 + outerRadius^2;
IMoIatO(1,1) = M*sumOfSquares/2;
IMoIatO(2,2) = (M/12)*(4*length^2 + 3*sumOfSquares);
IMoIatO(3,3) = IMoIatO(2,2);

% Transfer this mass moment of inertia matrix to the center of mass
IMoIatCoM = translateIMoIatPToCoM(M, CoM, IMoIatO);

% Compute the time derivative of the inertia matrix if the motor is still
% burning. The motor is still burning if r0 < rp < r1
if (innerRadiusatt < innerRadius || innerRadiusatt > outerRadius)
    error("Unphysical burning.")
elseif (innerRadius == innerRadiusatt || innerRadiusatt == outerRadius)
    dWdt = 0;
else
    dWdt = (outerRadius - innerRadius)/burnTime;
end

% Compute the time derivative of the inertia matrix
dIMoIatCoMdt = 6*innerRadiusatt*dWdt*[2, 0, 0; 0, 1, 0; 0, 0, 1];

% Provide the motor's (geometric) profile (x, outersurf, innersurf, density)
profile.x = [0; length];
profile.outersurf = [1; 1]*outerRadius;
profile.innersurf = [1; 1]*innerRadiusatt;
