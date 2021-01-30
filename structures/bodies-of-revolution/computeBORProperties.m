function [M, CoM, IMoIatCoM] = computeBORProperties(x, f, g, density)
% 
% Matt Werner (m.werner@vt.edu) - Jan 15, 2021
% 
% Compute the mass, center of mass, and mass moment of inertia of a (B)ody
% (O)f (R)evolution (BOR) described by the 2D functions f(x) and g(x) by
% numerical integration. The 3D shape generation is realized by rotating
% the two provided functions about the x-axis.
% 
% The outside surface of the body is defined by the f(x), and the
% inside surface of the body is defined by the g(x). The physical
% (realistic) body may be hollow (exhibiting an empty cavity) if 0 < g < f
% for all values x. The body is solid if g is equivalently 0.
% 
% The center of mass (CoM) and mass moment of inertia matrix (IMoIG) are
% computed with respect to the body's center of mass using a coordinate
% system exhibiting the the same orientation of the coordinate system
% defining f(x) and g(x) but simply shifted over to the location of the
% center of mass.
% 
% Note: The coordinate system defining f and g must be oriented such that
% y = f(x), y = g(x), and x is defined on the closed interval [0, L], where
% L > 0 is length of the body.
% 
%   y |
%     |
%     |
%     |          empty                                   ________  f(x)
%     |        space                           ,________/
%     |     (outer)                  ,________/            ,_____  g(x)
%     |                       ,_____/               ,_____/
%     |                 ,____/                ,____/
%     |            ,___/                 ,___/
%     |        ,__/       solid      ,__/           empty  
%     |    _,_/         space    _,_/             space
%     | _,/                   ,_/              (inner)
%     |/_____________________/__________|_______________________| x (revolve)
% z (.)                                CoM                      L
% 
%    Inputs:
% 
%                 x - Sampled points where the functions f(x) and g(x) are
%                     evaluated.
%                     Size: 1-by-n (vector)
%                     Units: m (meters)
% 
%                 f - Function representing the outer surface of the body
%                     of revolution.
%                     Size: 1-by-n (vector)
%                     Units: m (meters)
% 
%                 g - Function representing the inner surface of the body
%                     of revolution.
%                     Size: 1-by-n (vector)
%                     Units: m (meters)
% 
%            density - Material density function of the body of revolution
%                     at each sampled point x. If specified as a scalar,
%                     then the density is assumed uniform.
%                     Size: 1-by-n (vector) OR 1-by-1 (scalar)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%    Outputs:
% 
%                 M - The BOR's total mass.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               CoM - The distance of the BOR's center of mass from the
%                     coordinate origin.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%             IMoIG - The BOR's mass moment of inertia matrix taken with
%                     respect to the BOR's center of mass using a
%                     coordinate system whose orientation is the same as
%                     that used to defined the BOR's profile with f and g.
%                     Size: 3-by-3 (matrix)
%                     Units: kg*m2 (kilograms times squared meters)
% 

% Precompute
f2 = f.^2;
g2 = g.^2;
x2 = x.^2;

% Define the auxillary functions A (local area), B, and C
A = pi*(f2 - g2);
B = f2 + g2;
C = B/2 + 2*x2;

% Compute the mass per unit length
MPUL = density.*A;

% Compute the mass
M = trapz(x, MPUL);

% Compute the center of mass
CoM = [trapz(x, MPUL.*x)/M; 0; 0];

% Compute the mass moment of inertia matrix relative to the coordinate
% frame used to define f and g (same orientation and origin)
IMoIatOrigin = diag(trapz(x, MPUL .* [B, C, C] / 2));

% Compute the mass moment of inertia matrix relative to the coordinate
% frame whose orientation is the same as that of the frame used to define
% the f and g functions but simply shifted over to have its origin located
% at the BOR's center of mass. This transformation is achieved via the 3D
% parallel axis theorem
IMoIatCoM = IMoIatOrigin - M*(CoM'*CoM*eye(3) - CoM*CoM');