function Reff = computeEffectiveEarthRadius(geodeticLat, Req, f)
% 
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Calculate an effective (spherical) Earth's radius at a geodetic latitude
% with respect to the physical Earth's ellipsoid. The effective radius is
% the geometric distance between the center of the ellipsoid and the point
% on the ellipsoid's surface corresponding to the specified geodetic
% latitude (at any arbitrary longitude). The effective radius in terms of
% the geodetic latitude is written
%                     _________________________________
%                    / cos2(lat) + (1 - f)^4 sin2(lat)
%     R    = R   \  / ---------------------------------,
%      eff    eq  \/   cos2(lat) + (1 - f)^6 sin2(lat)
% 
% where sin2(lat) = sin(lat)^2 and cos2(lat) = cos(lat)^2 with lat being
% the geodetic latitude.
% 
% 
%    Inputs:
% 
%       geodeticLat - Geodetic latitude of the current position with
%                     respect to the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%               Req - Equatorial (semimajor) radius of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
% 
%                 f - Flattening of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 

% Precalculate (1 - f)^4 and (1 - f)^6 using the facts that (-1)^2 = 1 and
% (x^a)^b = x^(a*b)
fminus1ToThe4 = (1 - f)^4;
fminus1ToThe6 = fminus1ToThe4^1.5;
% Precalculate cos(lat)^2 and sin(lat)^2
cos2lat = cos(geodeticLat)^2;
sin2lat = sin(geodeticLat)^2;

% Calculate the effective radius
Reff = Req * sqrt((cos2lat + fminus1ToThe4*sin2lat) / (cos2lat + fminus1ToThe6*sin2lat));