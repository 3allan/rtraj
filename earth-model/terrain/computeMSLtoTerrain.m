function H = computeMSLtoTerrain(longitude, geocentricLatitude, HC, HS)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the height (orthometric) from mean sea level (MSL) defined by the
% WGS84 ellipsoid to the real/actual ground (terrain) above it.
% 
%    Inputs: 
% 
%          longitude - Longitude of the current position on (or
%                      above/below) the WGS84 ellipsoid, which is expressed
%                      in spherical coordinates relative to the Earth-
%                      Centered-Fixed (ECF) frame, wherein the +x
%                      axis is defined to be located at 0 degrees longitude
%                      and defines a meridian arc passing directly through
%                      Greenwich, England and the +z axis passes through
%                      the North pole.
%                      Size: 1-by-1 (scalar)
%                      Units: - (radians)
% 
% geocentricLatitude - Latitude of the current position on (or above/below)
%                      the WGS84 ellipsoid, but given in terms of the 
%                      geocentric (not geodetic) latitude, which is 
%                      expressed in spherical coordinates relative to the
%                      Earth-Centered-Fixed (ECF) frame, wherein the +x
%                      axis is defined to be located at 0 degrees longitude
%                      and defines a meridian arc passing directly through
%                      Greenwich, England and the +z axis passes through
%                      the North pole. Note that this angle is not the
%                      geocentric colatitude.
%                      Size: 1-by-1 (scalar)
%                      Units: - (radians)
% 
%                 HC - Optional(!) Cosine-series harmonic coefficients for
%                      the elevation of local topography ("the ground")
%                      relative to the EGM08 geoid, which approximates the
%                      notion of mean sea level (MSL) at all points on the
%                      globe (even over land). This kind of height is
%                      called the orthometric height. If HC is specified,
%                      then HS must also be specified.
%                      Size: 2401336-by-1 (vector)
%                      Units: m (meters)
% 
%                 HS - Optional(!) Sine-series harmonic coefficients for
%                      the elevation of local topography ("the ground")
%                      relative to the EGM08 geoid, which approximates the
%                      notion of mean sea level (MSL) at all points on the
%                      globe (even over land). This kind of height is
%                      called the orthometric height. If HS is specified,
%                      then HC must also be specified and the orthometric
%                      height will be calculated using the harmonic series
%                      approach.
%                      Size: 2401336-by-1 (vector)
%                      Units: m (meters)
% 
%    Outputs:
% 
%                  H - Orthometric height of the EGM08 geoid relative to
%                      the WGS84 ellipsoid. Orthometric height in this
%                      context is equivalent to the concept of ellipsoidal
%                      height at a given longitude and geodetic latitude.
%                      Size: 1-by-1 (scalar)
%                      Units: m (meters)
% 

% Check coefficients
checkCoefficients(HC, HS)

% Precalculate the argument of the associated Legendre polynomial
sinGeocentricLatitude = sin(geocentricLatitude);
cosLon = cos(longitude);
sinLon = sin(longitude);

% Obtain information about the degree of the system (should be 2190 using
% the Coeff_Height_and_Depth_to2190_DTM2006.0 file provided by the NGA as
% part of EGM2008.
L = numel(HC);
if (L == 2401336)
    % Specific number for complete degree and order of 2190
    N = 2190;
else
    % Maximum degree of complete degree and order system starting at (0,0).
    N = floor(0.5*(-3 + sqrt(1 + 8*L)));
end

% Compute the series
H = 0; % Orthometric height
nmidx = 1;
for nn = 0:N
    % Precompute all necessary Legendre polynomials for this degree
    PnmAll = computefnLegendre(nn, sinGeocentricLatitude);
    for mm = 0:nn
        % Increase efficiency by utilizing recurrence relations on
        % trigonometric functions
        if (mm > 0)
            cosmLon = cosmLonNext;
            sinmLon = sinmLonNext;
        else
            cosmLon = 1;
            sinmLon = 0;
        end
        % Obtain these coefficients of degree n and order m amidst the tall
        % and well-ordered stack
        HCnm = HC(nmidx);
        HSnm = HS(nmidx);
        Pnm = PnmAll(1+mm);
        H = H + (HCnm*cosmLon + HSnm*sinmLon)*Pnm;
        % Recurrence relations
        cosmLonNext = cosLon*cosmLon - sinLon*sinmLon;
        sinmLonNext = sinLon*cosmLon + cosLon*sinmLon;
        % Increase count
        nmidx = nmidx + 1;
    end
end