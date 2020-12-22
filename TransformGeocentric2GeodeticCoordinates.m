function [longitude, geodeticLat, GPSh] = TransformGeocentric2GeodeticCoordinates(x, y, z, Req, e)
% 
% Matt Werner (m.werner@vt.edu) - Dec 22, 2020
% 
%    Inputs:
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
% 

% Compute the longitude exactly
longitude = atan2(y, x);

% Calculate some repetetive and constant quantities
p = sqrt(x^2 + y^2);
e2 = e^2;

% Compute the latitude and height above the ellipsoid
if (p < Req*1.5679e-06)
    % Obtain results approximately exactly (positioned close to directly
    % over either the North or South pole). 
    % Note: 100 = 63781370*1.5679e-06
    Rpo = Req*sqrt(1 - e^2);
    GPSh = abs(z) - Rpo;
    geodeticLat = sign(z)*pi/2;
    return
else
    % Iteratively obtain the results (not positioned over the poles)
    geodeticLatTmp = atan(z / ((1 - e^2)*p));
    while true
        RnTmp = computePrimeVerticalRadius(Req, e, geodeticLatTmp, "geodetic");
        GPShTmp = p*sec(geodeticLatTmp) - RnTmp;
        geodeticLatOld = geodeticLatTmp;
        geodeticLatTmp = atan(z / ((1 - e2*RnTmp / (RnTmp + GPShTmp))*p));
        if (norm(geodeticLatTmp - geodeticLatOld) < 1e-6)
            % Assign converged (geodetic) latitude and (ellipsoidal) height
            geodeticLat = geodeticLatTmp;
            GPSh = GPShTmp;
            return
        end
    end
end