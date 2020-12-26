function [x, y, z] = TransformGeodetic2GeocentricCoordinates(GPSh, ...
                        geodeticLat, longitude, Req, e)
% 
% Matt Werner (m.werner@vt.edu) - Dec 22, 2020
% 
%    Inputs:
% 
%              GPSh - Height above the surface of the ellipsoid at a
%                     particular geodetic latitude and longitude. Here,
%                     "above" indicates that the height be measured normal
%                     to the surface of the ellipsoid starting from the
%                     point on the ellipsoid's surface described by the
%                     angular coordinates of geodetic latitude and
%                     longitude.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%       geodeticLat - Geodetic latitude relative to the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: rad (radians)
% 
%         longitude - Longitude relative to the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: rad (radians)
% 
%               Req - Equatorial radius (semimajor axis) of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 e - Eccentricity of the ellipsoid defined by the standard
%                     relationship with the flattening wherein
%                                       e2 = 2f - f2,
%                     where e is taken along the positive branch of the
%                     square root.
%                     Size: 1-by-1 (scalar)
%                     Units: - (N/A)
% 
%    Outputs:
% 
%                 x - Geocentric position measured from the center of Earth
%                     along the equatorial plane in the direction of
%                     precisely where 0 degrees longitude is defined
%                     (Greenwich, England).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 y - Geocentric position measured from the center of Earth
%                     along the equatorial plane in the direction of
%                     precisely where 90 degrees longitude is defined
%                     (in the Bay of Bengal).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 z - Geocentric position measured from the center of Earth
%                     normal to the equatorial plane in the direction of
%                     the North pole.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 

% Check that the specified latitude and longitude are in the correct range
checkxInInterval(geodeticLat, [-pi/2, pi/2])
checkxInInterval(longitude, [-pi, pi])

% Compute the prime vertical radius of the ellipsoid (having an associated
% equatorial (semimajor) radius Req and eccentricity e) at the specified 
% geodetic latitude.
Rn = computePrimeVerticalRadius(Req, e, geodeticLat, "geodetic");

% Precompute repeated quantities
Rn_plus_GPSh = Rn + GPSh;
cosLat = cos(geodeticLat);

% Define the transformation from geodetic coordinates to geocentric
% coordinates
x = Rn_plus_GPSh*cosLat*cos(longitude);
y = Rn_plus_GPSh*cosLat*sin(longitude);
z = ((1 - e^2)*Rn + GPSh)*sin(geodeticLat);