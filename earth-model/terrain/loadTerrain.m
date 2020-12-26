function [longitude, geodeticLatitude, WGS84ToGeoid, GeoidToTerrain, ...
          WGS84ToTerrain] = loadTerrain(angularUnits)
% 
% Matt Werner (m.werner@vt.edu) - Dec 2, 2020
% 
% Load terrain models according to the specified units for the longitude
% and geodetic latitude angles.
% 
%    Inputs:
% 
%      angularUnits - Optional(!) Specifies which units are desired for the
%                     representations of longitude and geodetic 
%                     latitude determining how far the geoid (defined by
%                     EGM2008) and local topography are from the WGS84
%                     ellipsoid as well as how far the geoid is from the
%                     ellipsoid. With the geoid representing mean sea level
%                     (MSL) globally and ellipsoidal height determined
%                     globally by GPS, the orthometric height (height above
%                     MSL) may be determined.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                     "radians" - Provides the longitude and geodetic
%                                 latitude in units of radians (default).
% 
%                     "degrees" - Provides the longitude and geodetic
%                                 latitude in units of degrees.
% 
%    Outputs:
% 
%         longitude - Longitude associated with Earth-Centered-Fixed (ECF)
%                     coordinates. The longitude is defined by ellipsoidal
%                     coordinates (which, coincidentally, shares the same
%                     definition with spherical coordinates) being
%                     identically zero for its meridian arc to pass
%                     directly through Greenwich, England and is
%                     represented by negative values (as (-pi (or -180), 0)
%                     radians (or degrees)) for all locations west of 
%                     Greenwich and positive values (as (0, +pi (or +180))
%                     radians (or degrees)) for all locations east of
%                     Greenwich. Thus, the longitude is represented in the
%                     interval (-pi, pi) or (-180, 180) degrees depending
%                     on the desired units.
%                     Size: 1-by-N (vector)*
%                     Units: - (radians or degrees)
% 
%  geodeticLatitude - Geodetic (not geocentric) latitude associated with
%                     Earth-Centered-Fixed (ECF) coordinates. Geodetic
%                     latitude is defined by ellipsoidal coordinates
%                     with zero latitude corresponding to the line of
%                     points along the equator (the great circle) and is
%                     represented by negative values (as (-pi/2 (or -90), 0)
%                     radians (or degrees)) for all locations south of the
%                     equatorial plane and positive values (as (0, +pi/2)
%                     (or +90)) radians (or degrees)) for all locations
%                     north of the equatorial plane. Thus, the geodetic
%                     latitude is represented in the interval (-pi/2, pi/2)
%                     or (-90, 90) degrees depending on the desired units.
%                     Again, geodetic latitude for an ellipsoid is NOT 
%                     defined by the usual spherical coordinate 
%                     representation of latitude/colatitude.
%                     Size: N-by-1 (vector)*
%                     Units: - (radians or degrees)
% 
%      WGS84ToGeoid - Height (ellipsoidal) of the EGM2008 geoid relative
%                     to the WGS84 ellipsoid worldwide determined by
%                     reading the 2.5'-by-2.5' resolution big-endian map
%                     computed using EGM2008 tide-free provided by the
%                     National Geospace Intelligence Agency (NGA) Office of
%                     Geomatics. File is available for download at
%                     https://earth-info.nga.mil/GandG/wgs84/gravitymod/...
%                       ...egm2008/egm08_wgs84.html
%                     Size: N-by-M (matrix)*
%                     Units: m (meters)
% 
%    GeoidToTerrain - Height (orthometric) of local topography ("the
%                     ground") relative to the EGM2008 geoid, which
%                     represents mean sea level (MSL) everywhere, even over
%                     land. These values were determined with spherical
%                     harmonic analysis using the elevation coefficients
%                     HCnm and HSnm provided by the National Geospace
%                     Intelligence Agency (NGA) Office of Geomatics. The
%                     resolution of the grid points was chosen to mimic the
%                     provided file used to create WGS84ToGeoid. The
%                     resolution is 2.5'-by-2.5'. File containing the
%                     coefficients is avilable for download at
%                     https://earth-info.nga.mil/GandG/wgs84/gravitymod/...
%                       ...egm2008/first_release.html
%                     Size: N-by-M (matrix)*
%                     Units: m (meters)
% 
%    WGS84ToTerrain - Height (ellipsoidal) of local topography ("the
%                     ground") relative to the WGS84 ellipsoid. These
%                     values are obtained according to
%                        WGS84ToGround = WGS84ToGeoid + GeoidToGround,
%                     provided that the grid points for each of the two
%                     elements on the right-hand-side are identically the
%                     same. WGS84ToGround values at a given longitude and
%                     geodetic latitude will, in general, NOT match what
%                     GPS satellites report the ellipsoidal height to be at
%                     that same location (latitude, geodetic latitude).
%                     This fact is due in part because of the 2.5'-by-2.5'
%                     resolution of the two grids, which is the lesser
%                     option of the two provided by NGA (the other is
%                     1'-by-1' resolution), and also due to interpolation
%                     errors.
%                     Size: N-by-M (matrix)*
%                     Units: m (meters)
% 
%                   * With 2.5'-by-2.5' grid resolution, one has that
%                     N = 8640 and M = 4321
%                     (8640 = 360/(2.5/60) and 4321 = 180/(2.5/60) + 1)
% 

% Display
disp("Loading terrain...")

% Check that the input is of type string/char
checkInput(angularUnits)

% Lowercase
angularUnits = lower(angularUnits);

if (nargin == 1)
    switch angularUnits
        case "radians"
            load("mats/EGM/wgs84-to-egm08geoid_undulation=meters_LatLon=(egm08geoid-to-ground_LatLon=radians).mat", ...
                "LongitudeW2GwithInterpolatedG2Tcoords", "geodeticLatitudeW2GwithInterpolatedG2Tcoords", ...
                "W2GwithInterpolatedG2Tcoords");
            load("mats/EGM/egm08geoid-to-ground_orthometricHeight=meters_LatLon=radians.mat", ...
                "G2T")
            load("mats/EGM/wgs84-to-ground_ellipsoidalHeight=meters_LatLon=radians.mat", ...
                "W2T")
        case "degrees"
            load("mats/EGM/wgs84-to-egm08geoid_undulation=meters_LatLon=(egm08geoid-to-ground_LatLon=degrees).mat", ...
                "LongitudeW2GwithInterpolatedG2Tcoords", "geodeticLatitudeW2GwithInterpolatedG2Tcoords", ...
                "W2GwithInterpolatedG2Tcoords")
            load("mats/EGM/egm08geoid-to-ground_orthometricHeight=meters_LatLon=degrees.mat", ...
                "G2T")
            load("mats/EGM/wgs84-to-ground_ellipsoidalHeight=meters_LatLon=degrees.mat", ...
                "W2T")
        otherwise
            error("Invalid angle units.")
    end
else
    % Load default units (radians)
    load("mats/EGM/wgs84-to-egm08geoid_undulation=meters_LatLon=(egm08geoid-to-ground_LatLon=radians).mat", ...
                "LongitudeW2GwithInterpolatedG2Tcoords", "geodeticLatitudeW2GwithInterpolatedG2Tcoords", ...
                "W2GwithInterpolatedG2Tcoords");
    load("mats/EGM/egm08geoid-to-ground_orthometricHeight=meters_LatLon=radians.mat", ...
        "G2T")
    load("mats/EGM/wgs84-to-ground_ellipsoidalHeight=meters_LatLon=radians.mat", ...
        "W2T")
end

% Distribute results
longitude = LongitudeW2GwithInterpolatedG2Tcoords;
geodeticLatitude = geodeticLatitudeW2GwithInterpolatedG2Tcoords;
WGS84ToGeoid = W2GwithInterpolatedG2Tcoords;
GeoidToTerrain = G2T;
WGS84ToTerrain = W2T;

% Display
disp("Loaded terrain")