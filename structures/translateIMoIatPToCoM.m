function IMoIatCoM = translateIMoIatPToCoM(M, CoMToP, IMoIatP)
% 
% Matt Werner (m.werner@vt.edu) - Jan 31, 2021
% 
% Shift the expression of the mass moment of inertia matrix of a rigid body
% from an arbitrary point P to the body's center of mass.
% 
%    Inputs:
% 
%                 M - Rigid body mass whose inertia matrix is being shifted
%                     over to its center of mass.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%            CoMToP - Position of the body's center of mass relative to
%                     some other point P.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%           IMoIatP - Mass moment of inertia matrix of the rigid body
%                     relative to some other point P.
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%    Outputs:
% 
%         IMoIatCoM - Mass moment of inertia matrix of the rigid body
%                     relative to the rigid body's center of mass. The
%                     orientation of the frame expressing this matrix is
%                     exactly the same as the orientation expressing the
%                     center of mass relative to the point P (CoM) and the
%                     mass moment of inertia matrix relative to the point P
%                     (MoI).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 

% Transfer the mass moment of inertia matrix about the point P to about the
% rigid body's center of mass
IMoIatCoM = IMoIatP - M*(CoMToP'*CoMToP*eye(3) - CoMToP*CoMToP');