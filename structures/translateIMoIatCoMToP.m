function IMoIatP = translateIMoIatCoMToP(M, PToCoM, IMoIatCoM)
% 
% Matt Werner (m.werner@vt.edu) - Jan 31, 2021
% 
% Shift the expression of the mass moment of inertia matrix of a rigid body
% from the body's center of mass to an arbitrary point P.
% 
%    Inputs:
% 
%                 M - Rigid body mass whose inertia matrix is being shifted
%                     over to the arbitrary point P.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%            PToCoM - Position of the arbitrary point P relative to the
%                     rigid body's center of mass.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%         IMoIatCoM - Mass moment of inertia matrix of the rigid body
%                     relative to the body's center of mass.
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 
%    Outputs:
% 
%           IMoIatP - Mass moment of inertia matrix of the rigid body
%                     relative to the rigid body's center of mass. The
%                     orientation of the frame expressing this matrix is
%                     exactly the same as the orientation expressing the
%                     center of mass relative to the point P (CoM) and the
%                     mass moment of inertia matrix relative to the point P
%                     (MoI).
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 

% Transfer the mass moment of inertia matrix about the rigid body's center
% of mass to about the point P
IMoIatP = IMoIatCoM + M*(PToCoM'*PToCoM*eye(3) - PToCoM*PToCoM');