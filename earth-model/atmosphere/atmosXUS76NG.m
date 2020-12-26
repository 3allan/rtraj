function [d, T, p, c, mu] = atmosXUS76NG(geomH, geodeticLatitude, Req, f)
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
%                     surface. The geopotential altitude is calculated
%                     according to the spherical Earth of effective radius
%                     Reff.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
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
%                 p - Atmospheric pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%         varargout - Additional outputs specified.
%                     Size: ?
%                     Units: ?
% 

% Define the Extended 1976 U.S. Standard Atmosphere, but using an effective
% Earth radius and normal gravitational acceleration
Reff = computeEffectiveEarthRadius(geodeticLatitude, Req, f);
geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH);
gSL = computeNormalGravity(geodeticLatitude, 0, Req, f);

% Access the Extended 1976 U.S. Standard Atmosphere using this slightly
% different definition of geopotential altitude and gravitational
% acceleration at sea level, which are both dependent upon the geodetic
% latitude of position.
[d, T, p, c, mu] = atmosXUS76(geopH, gSL, Reff);