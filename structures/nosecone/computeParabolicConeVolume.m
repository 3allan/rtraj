function volume = computeParabolicConeVolume(radius, length, thickness, k)
% 
% Matt Werner (m.werner@vt.edu) - Jan 14, 2021
% 
% Reference the mass of a parabolic cone. The cone may be completely solid
% or completely hollow, consisting only of a shell of material on its
% outside. Determining the cone's hollowness is determined by the
% thickness, where 'length' thickness corresponds with a solid cone and 0 
% thickness corresponds with the shell.
% 
% Parabolic cones do include the linear variety that have no curvature to
% their sides, indicated by a cone constant k = 0. Otherwise, the cone has
% curvature due to a cone constant 0 < k < 1, where k = 1 corresponds to a
% cone that has vanishing slope at its base (smoothly transitions to the
% body).
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
%    Outputs:
% 
%            volume - The space occupied by the solid structure defining
%                     the cone. If the structure is constructed such that
%                     its material density is constant throughout, then
%                     this volume multiplied with the material density
%                     provides the cone's mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m3 (cubic meters)
% 

% Checks
x = [radius; length; thickness; k];
I = [eps, inf; eps, inf; 0, length; 0, 1];
checkxInInterval(x, I)

% Define nondimensional parameters
l = thickness/length;

% Precompute repeated factors
l2 = l^2;
l3 = l2*l;
l4 = l3*l;
k2 = k^2;

% Nondimensional factor accounting for the cone's hollowness and curvature,
% dependent upon the thickness and curvature (k) parameters
num = l*(20*(l2 - 3*l + 3) + 15*(l3 - 4*l2 + 6*l - 4)*k + 3*(l4 - 5*l3 + 10*l2 - 10*l + 5)*k2);
den = 5*(2 - k)^2;
correctionFactor = num/den;

% Volume of a solid circular cone
volumeCircularCone = pi*radius^2*length/3;

% Calculate the volume of a parabolic cone as a perturbation away from the
% volume of a solid circular cone of equivalent radius and length
volume = volumeCircularCone*correctionFactor;