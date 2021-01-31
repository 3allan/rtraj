function [newM, newCoM, newIMoI] = attachComponents(M, CoM, IMoI, ...
                                    newMs, newCoMs, newIMoIs)
% 
% Matt Werner (m.werner@vt.edu) - Jan 30, 2021
% 
% Calculate the new position of a rigid body's center of mass due to adding
% new rigid-body components whose masses and centers of mass are provided
% and combined with the pre-existing body's mass and center of mass.
% 
%    Inputs:
% 
%                 M - Rigid body's current mass before adding any
%                     additional components.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%               CoM - Rigid body's current center of mass before adding any
%                     additional components. The frame defining this center
%                     of mass and the mass moment of inertia matrix must be
%                     identical.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%              IMoI - Rigid body's mass moment of inertia matrix before
%                     adding any additional components. The frame defining
%                     this matrix and the center of mass must be identical.
% 
%             newMs - Masses of the new components to be added to the
%                     rigid body.
%                     Size: 1-by-n (vector)
%                     Units: kg (kilograms)
% 
%           newCoMs - Positions of the new components' centers of mass
%                     relative to the same coordinate frame defining the
%                     rigid body's center of mass before adding any
%                     additional components.
%                     Size: 3-by-n (matrix)
%                     Units: kg (kilograms)
% 
%           newMoIs - Expressions of the new components' mass moments of
%                     inertia relative to the same coordinate frame
%                     (including the origin) defining the rigid body's mass
%                     moment of inertia matrix before adding any additional
%                     components.
%                     Size: 3-by-3*n (matrix)
%                     Units: kg (kilograms)
% 
%    Outputs:
% 
%            newCoM - New position of the rigid body's center of mass after
%                     securing the additional components to the original
%                     body.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 

% Leave dimensional checks to Matlab

% Get the number of new components to add to the rigid body
numComponents = size(newMs, 2);

% Compute the new mass
newM = M + sum(newMs, 2);

% Compute the new center of mass
newCoM = (M*CoM + sum(newMs.*newCoMs, 2)) / newM;

% Compute the new mass moment of inertia matrix
for component = 1:numComponents
    idx = 1 + 3*(component - 1);
    IMoIc = newIMoIs(:, idx:2+idx);
    newIMoI = IMoI + IMoIc;
end