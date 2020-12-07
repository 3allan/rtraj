function [d, T, varargout] = atmosXUS76NG(geomH, nreqout, ...
                                          geodeticLatitude, Req, f)
% 
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Calculate the Jet, 1976 U.S. Standard Atmosphere, and Extended 1976 U.S.
% Standard Atmosphere models given the geometric (physical) height above
% mean sea level (MSL) as well as the Extended 1976 U.S. Standard
% Atmosphere model but using normal gravity and an effective Earth radius
% to account for the fact that Earth's gravity varies like the ellipsoid.
% 
%    Inputs:
% 
%             geomH - Geometric (physical) height above mean sea level
%                     (MSL). Given only geomH, the Jet/1976USSA atmospheric
%                     properties are calculated according to how many
%                     outputs are requested (nreqout).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%           nreqout - Optional(!) Number of requested outputs. By default,
%                     XUS76NG outputs only density (d) and environmental
%                     temperature (T) if nreqout is omitted or 2, but will
%                     also provide the speed of sound (c) and dynamic
%                     viscosity (mu) if 3 or 4, respectively.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 
%  geodeticLatitude - Optional(!) Geodetic latitude of a position on or
%                     above the Earth's surface as referenced by the
%                     ellipsoid. Specifying a geodetic latitude will
%                     choose a different reference radius and sea-level
%                     value for the gravitational acceleration than the
%                     defined values for the 1976 U.S. Standard Atmosphere.
%                     The reference radius will be chosen to be the radius
%                     from the ellipsoid's center to the point on the
%                     ellipsoid's surface (height 0) at the specified
%                     geodetic latitude. Correspondingly, the gravitational
%                     acceleration experienced at sea level will be chosen
%                     to be value of normal gravity on the ellipsoid's
%                     surface. Correspondingly, the geopotential altitude 
% 
%    Outputs:
% 
%                 d - Atmospheric density.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%                 T - Environmental (kinetic) temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%         varargout - Additional outputs specified.
%                     Size: ?
%                     Units: ?
% 

% Define the Extended 1976 U.S. Standard Atmosphere, but using an effective
% Earth radius and normal gravitational acceleration
Reff = computeEffectiveEarthRadius(geodeticLatitude, Req, f);
g0 = computeNormalGravity(geodeticLatitude, 0, Req, f);
% calculate geopotential height
% call atmosUS76
