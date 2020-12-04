function [GMsun, GMmoon] = getGMsOfSunAndMoon
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Obtain fixed values for the gravitational parameters (GM) of the sun and
% moon. The gravitational parameter of each body is the standard product of
% Newton's gravitational constant (G) and the mass of the body (M). Note
% that neither of G nor M are known very well independently due to how weak
% gravitational fields are and the inability to measure any celestial
% body's mass directly. Their product, however, constituting the
% gravitational parameter GM is known to much greater precision because
% they always appear together in the dynamics of spacecraft and other
% celestial bodies.
% 
%    Inputs:
%                   -
% 
%    Outputs:
% 
%             GMsun - Gravitational parameter (GM) of the sun.
%                     Size: 1-by-1 (scalar)
%                     Units: m3/s2 (cubic meters per squared seconds)
% 
%            GMmoon - Gravitational parameter (GM) of the moon.
%                     Size: 1-by-1 (scalar)
%                     Units: m3/s2 (cubic meters per squared seconds)
% 

GMsun = 1.327124400189e20; % Sun
GMmoon = 4.90486959e12; % Moon