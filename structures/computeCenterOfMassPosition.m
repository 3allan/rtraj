function Gpos = computeCenterOfMassPosition(Gposs, m, M)
% 
% Matt Werner (m.werner@vt.edu) - Dec 27, 2020
% 
% Calculate the center of mass of a system of particles. Such a system of
% particles in this context is taken to mean the system consisting of
% masses from individual rigid-body components centered at their own
% individual centers of mass relative to some reference coordinate system.
% Usually, this reference coordinate system is taken to be some sort of
% body-fixed frame stationary to the rigid-body.
% 
%    Inputs:
% 
%             Gposs - Locations of the centers of mass for each individual
%                     rigid-body component composing the whole rigid body.
%                     Size: 3-by-n (matrix)
%                     Units: m (meters)
% 
%                 m - Masses of the individual rigid-body components
%                     composing the whole rigid body.
%                     Size: 1-by-n (vector)
%                     Units: kg (kilograms)
% 
%                 M - Total mass of the whole rigid body.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 

% No checks

% Compute the location of the center of mass
Gpos = Gposs*m(:)/M;