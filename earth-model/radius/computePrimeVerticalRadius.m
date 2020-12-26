function Rn = computePrimeVerticalRadius(Req, e, latitude, geotype)
% 
% Matt Werner (m.werner@.vt.edu) - Dec 2 2020
% 
% Computes the prime vertical radius at a geodetic latitude for the
% ellipsoid model associated with the eccentricity e. The prime vertical
% radius is obtained according to the standard formula
%                 a / Rn = sqrt(1 - e2 sin2(x)),
% where Rn is the prime vertical radius, a is the equatorial (semimajor)
% radius of the ellipsoid associated with the eccentricity e, and x is the  
% geodetic latitude at which to evaluate the prime vertical radius.
% 
%    Inputs:
% 
%               Req - Equatorial radius (semimajor axis) of the ellipsoid
%                     model associated with the given eccentricity e.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 e - Eccentricity of the ellipsoid model associated with
%                     equatorial radius (semimajor axis) Req.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%          latitude - Latitude at which to evaluate the prime vertical 
%                     radius, Rn. The latitude varies in the interval 
%                     [-pi/2, pi/2] radians such that a latitude of 0 
%                     corresponds with a line of points (a great circle) 
%                     around the equator and a geodetic latitude of +pi/2
%                     radians corresponds with a single point on the North 
%                     pole of the ellipsoid (-pi/2 radians on the South
%                     pole  of the ellipsoid).
%                     Size: n-by-1 (vector)
%                     Units: - (radians)
% 
%           geotype - Specifies the type of latitude being supplied by
%                     `latitude' such that the correct expression for the
%                     prime vertical radius, Rn, is returned.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                     "geocentric" - Indicates that the latitude be taken
%                                    in the sense of spherical coordinates
%                                    with respect to Earth-Centered-Fixed
%                                    (ECF) coordinates, where the x-axis
%                                    passes directly through Greenwich,
%                                    England and the z-axis points through
%                                    the North pole. Note that the
%                                    geocentric latitude is not the same
%                                    as geocentric colatitude. Geocentric 
%                                    latitude measures angles relative to 
%                                    the equatorial plane (in which
%                                    latitude is defined to be zero)
%                                    whereas geocentric colatitude measures
%                                    angles relative to the axis passing
%                                    through the North pole (through which
%                                    the colatitude is defined to be zero).
% 
%                       "geodetic" - Indicates that the latitude be taken
%                                    in the sense of ellipsoidal
%                                    coordinates with respect to Earth-
%                                    Centered-Fixed (ECF) coordinates,
%                                    where the x-axis passes directly
%                                    through Greenwich, England and the
%                                    z-axis points through the North pole.
%                                    Geodetic coordinates (latitude,
%                                    longitude, and height) are those
%                                    coordinates reported by all GPS
%                                    satellites for any location on Earth.
% 
%    Outputs:
% 
%                Rn - Prime vertical radius at the specified geodetic
%                     latitude(s) of the ellipsoid model with equatorial
%                     radius (semimajor axis) a and eccentricity e.
%                     Size: n-by-1 (scalar)
%                     Units: m (meters)
% 

% Check that values are valid
checkxInInterval(e, [0, 1])
checkxInInterval(latitude, [-pi/2, pi/2])
% Check the string input
checkInput(geotype)
% Lowercase
geotype = lower(geotype);

% Compute the prime vertical radius at the specified latitude according to
% geotype.
switch geotype
    case "geodetic"
        Rn = Req ./ sqrt(1 - (e*sin(latitude))^2);
    case "geocentric"
        cos2lat = cos(latitude).^2;
        sin2lat = sin(latitude).^2;
        Rn = Req * sqrt((cos2lat + (1 - e^2)^2 * sin2lat) ./ (cos2lat + (1 - e^2)^3 * sin2lat));
    otherwise
        error("Please specify a valid geotype.")
end