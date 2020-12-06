function [g0, r0, p0, T0, Y] = getUS76BaseValues
% 
% Matt Werner (m.werner@vt.edu) - Dec 5, 2020
% 
% Get defining parameters for dry air at sea level in accordance with the
% U.S. Standard Atmosphere, 1976 atmospheric model.
% Reference: NASA-TM-X-74335 (Document ID 19770009539)
%  https://ntrs.nasa.gov/citations/19770009539
% 
%    Inputs:
%                   -
% 
%    Outputs:
% 
%                g0 - The sea-level value of the acceleration of gravity.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per second squared)
% 
%                r0 - The effective Earth's radius at the latitude at which
%                     the sea-level value of the acceleration of gravity,
%                     g0, is "precisely" 9.80665 m/s2.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters per second squared)
% 
%                p0 - The sea-level value of the atmospheric pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per second squared)
% 
%                T0 - The sea-level value of the air temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%                 Y - Ratio of specific heats for dry air.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)

% Provide defining parameters for the 1976 U.S. Standard Atmosphere
g0 = 9.80665;
r0 = 6356766;
p0 = 101325;
T0 = 288.15;
Y = 1.40;
