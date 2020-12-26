function [Y0, Y] = computeNormalGravity(geodeticLatitude, h, Req, f)
% 
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Compute the gravitational acceleration with respect to the ellipsoid
% defined by the equatorial (semimajor) radius Req and flattening f at the
% height above the ellipsoid h corresponding at a geodetic latitude as a
% quadratic function of h. This approximation is decent for ellipsoidal
% heights h much less than the equatorial radius (h << Req); it begins to
% diverge at altitudes such that h/Req ~ 0.2. For low-Earth orbits (VLEOs),
% h/Req ~ 0.03, which is well within the region of convergence.
% For reference, GPS satellites orbit at an altitude of about 20,000 km, or
% h/Req ~ 3.175, whereas the ISS orbits at an altitude of about 410 km, or
% h/Req ~ 0.064. The direction of normal gravity is, of course, normal to
% the ellipsoid, so in an East-North-Vertical coordinate system situated at
% a height h corresponding with the geodetic latitude and any (arbitrary)
% longitude, the gravitational acceleration is [0, 0, -Y]'.
% 
%    Inputs:
% 
%  geodeticLatitude - Geodetic latitude of the current position
%                     corresponding with the height above the ellipsoid.
%                     Size: 1-by-1 (string)
%                     Units: - (radians)
% 
%                 h - Height above the ellipsoid corresponding with the
%                     geodetic latitude. If h = 0, then the gravitational
%                     acceleration returned is the normal gravity at sea
%                     level (on the surface of the ellipsoid). Otherwise, a
%                     second-order Taylor-series expansion of the normal
%                     gravity is returned, which is good for altitude even
%                     exceeding low-Earth orbit (LEO) a little bit.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               Req - Equatorial radius (semimajor axis) of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 f - Flattening of the ellipsoid.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%    Outputs:
% 
%                Y0 - Normal gravity to the ellipsoid at the specified
%                     geodetic latitude at sea level (on the ellipsoid's
%                     surface). The normal gravity in this context is 
%                     specified to be positive.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per squared second)
%
%                 Y - Normal gravity to the ellipsoid at the specified
%                     geodetic latitude and ellipsoidal height. The normal
%                     gravity in this context is specified to be positive.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per squared second)

% Precalculate sin^2(geodeticLatitude)
sin2lat = sin(geodeticLatitude)^2;

% Compute normal gravity at sea level
% Y0, the normal gravity at sea level, takes the form
%                       1 + 0.00193185138639*sin2lat
% Y0 = -9.7803267714 ------------------------------------
%                     (1 - 0.00669437999013*sin2lat)^0.5
% as per the National Geospatial-Intelligence Agency (NGA). The top two
% coefficients have been premultiplied in the following expression as to
% not waste computation power and the minus sign removed to ensure that the
% resulting normal gravity is positive to be consistent with the fact that
% the norm of the vector Yvec = [0, 0, -Y] is |Yvec| = Y > 0.
Y0 = (9.7803267714 + 0.0188941378326763*sin2lat) / sqrt(1 - 0.00669437999013*sin2lat);

% Check if an expansion is needed
if (nargin == 1 || h == 0)
    Y = Y0;
else
    % Precalculate a repeating term
    hOverReq = h/Req;
    % Check that a reasonable ellipsoidal height has been supplied
    if (hOverReq > 0.2)
        error("Normal gravity diverges from second-order approximation for requested altitude (h/Req ~ %1.2f)", hOverReq)
    end
    Y = Y0 * (1 - 2*(1 + f - 2*f*sin2lat)*hOverReq + 3*hOverReq^2);
end