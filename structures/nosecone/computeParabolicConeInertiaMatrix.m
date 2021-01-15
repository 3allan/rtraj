function IMoIG = computeParabolicConeInertiaMatrix(radius, length, ...
    thickness, k, mass, CoM)
% 
% Matt Werner (m.werner@vt.edu) - Jan 14, 2021
% 
% Calculate the mass moment of inertia matrix of a parabolic cone with
% respect to the cone's center of mass and coordinates directions such
% that, when placed at the cone's vertex, are defined:
%  - x: Points towards the center of mass. This direction is exactly
%       along the cone's line of symmetry pointing inwards from the cone's
%       vertex to its base.
%  - y: Perpendicular to the x-axis. Due to radial symmetry, the exact
%       direction of y is meaningless to specify.
%  - z: This direction completes the right-handed rule.
% 
% The moment of inertia for the parabolic has been precomputed using its
% exact integral definition and is reported with respect to the center of
% mass along the x,y,z axes defined precisely above.
% 
%    Inputs:
% 
%            radius - Radius of the cone at its base.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%            length - Length of the cone from its tip to its base. This
%                     quantity is otherwise known as the cone's height.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%         thickness - Distance between the physical cone's outside tip and
%                     inside tip, where the inside tip is realized by
%                     shifting a copy of the physical cone over by a
%                     distance equal to the thickness and taking its
%                     negative space. Because parabolic cones are convex
%                     shapes, this definition of thickness is well-defined
%                     for values between 0 and the indicated length.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 k - Cone constant determining the cone's curvature. A
%                     cone constant of zero results in a shape identifying
%                     precisely with a right circular cone having an
%                     infinitely sharp vertex. Increasing the cone constant
%                     upwards towards its upper limit of 1 increases the
%                     cone's curvature. In this case, the vertex is still
%                     sharp, but sides at least curve into it smoothly. A
%                     cone constant of 1 results in the cone having a
%                     vanishing derivative of its profile at its base such
%                     that the meeting between the side and base face
%                     occurs perpendicularly.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%              mass - Total mass of the cone.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%               CoM - Distance of the cone's center of mass from the cone's
%                     vertex.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%    Outputs:
% 
%             IMoIG - The cone's mass moment of inertia that describes the
%                     mass distribution relative to a reference point and
%                     coordinate system. The reference point here is taken
%                     to be the cone's center of mass and the coordinate
%                     system is defined above in the description. Due to
%                     radial symmetry, the inertia matrix is of the
%                     principal variety such that all products of inertia
%                     are identically zero.
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms * squared meters)
% 

% Checks
x = [length; thickness; k];
I = [eps, inf; 0, length; 0, 1];
checkxInInterval(x, I)

% Define nondimensional parameters
a = radius/length;
l = thickness/length;

% Precompute repeating factors
a2 = a^2;
l2 = l^2;
l3 = l2*l;
l4 = l3*l;
l5 = l4*l;
l6 = l5*l;
l7 = l6*l;
l8 = l7*l;
k2 = k^2;
k3 = k2*k;
k4 = k3*k;

%% Ixx (along the cone height direction)

% Inertia component of a shelled (completely hollow) right circular cone
% about the height direction
IxxShellCircularCone = mass*radius^2/2;

% Additional repeating terms
lm1 = l - 1;
lm15 = lm1^5;
lm16 = lm15*lm1;
lm17 = lm16*lm1;
lm18 = lm17*lm1;
lm19 = lm18*lm1;
km22 = (2 - k)^2;

% Nondimensional factor accounting for the cone's additional mass and
% curvature, dependent upon thickness (l) and curviture (k)
num = 1008*(1 + lm15) - 1680*(1 - lm16)*k + 1080*(1 + lm17)*k2 - 315*(1 - lm18)*k3 + 35*(1 + lm19)*k4;
den = 21 * (20*(l3 - 3*l2 + 3*l) + 15*(l4 - 4*l3 + 6*l2 - 4*l)*k + 3*(l5 - 5*l4 + 10*l3 - 10*l2 + 5*l)*k2) * km22;
correctionFactor = num/den;

% Calculate the Ixx component at the cone's center of mass along the
% coordinate frame defined above as a perturbation away from that of the
% right circular cone shell
Ixx = IxxShellCircularCone*correctionFactor;
clear IxxShellCircularCone correctionFactor

%% Iyy & Izz (perpendicular to the cone height direction)

% Inertia component of a shelled (completely hollow) right circular cone
% about the height direction but having its base radius the size as this
% cone's length (height)
IxxShellCircularCone = mass*length^2/2;

% Nondimensional factor accounting for the cone's additional mass and
% curvature, dependent upon thickness (l), wideness (a), and curviture (k)
num = +336*(2*(l4 - 10*l + 15) + 3*(l4 - 5*l3 + 10*l2 - 10*l + 5)*a2) ...
      -336*((-l5 + 2*l4 + 20*l2 - 65*lm1 + 1) + 5*(-l5 + 6*l4 - 15*l3 + 20*l2 - 15*l + 6)*a2)*k ...
      + 24*((2*l6 - 14*l5 + 7*l4 - 70*l3 + 490*l2 - 952*l + 749) + 45*(l6 - 7*l5 + 21*l4 - 35*l3 + 35*l2 - 21*l + 7)*a2)*k2 ...
      +  3*(-4*(4*l6 - 7*l5 - 140*l3 + 560*l2 - 819*l + 532) + 105*(l7 - 8*l6 + 28*l5 - 56*l4 + 70*l3 - 56*l2 + 28*l - 8)*a2)*k3 ...
      +    (12*(l6 - 35*l3 + 105*l2 - 126*l + 70) + 35*(l8 - 9*l7 + 36*l6 - 84*l5 + 126*l4 - 126*l3 + 84*l2 - 36*l + 9)*a2)*k4;
den = 42 * (20*(l2 - 3*l + 3) + 15*(l3 - 4*l2 + 6*l - 4)*k + 3*(l4 - 5*l3 + 10*l2 - 10*l + 5)*k2) * km22;
correctionFactor = num/den - 2*(CoM/length)^2; % Includes 3D parallel axis theorem using scalar CoM

% Calculate the Iyy & Izz components at the cone's center of mass along the
% coordinate frame defined above as a perturbation away from that of the
% right circular cone shell
Iyy = IxxShellCircularCone*correctionFactor;
Izz = Iyy;

IMoIG = [Ixx,   0,    0
          0,   Iyy,   0
          0,    0,   Izz];
