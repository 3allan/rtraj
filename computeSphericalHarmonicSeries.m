function SHS = computeSphericalHarmonicSeries (N, M, C, S, ...
                longitude, latitude, radius, A, B)
% 
% Matt Werner (m.werner@vt.edu) - Dec 3, 2020
% 
% Compute the spherically harmonic series obtained as the solution to
% Laplace's equation in spherical coordinates as a function of at least
% longitude and latitude. The series solution corresponding to the
% Laplace's equation
%          2        2        2
%         d  u     d  u     d  u 
%   0 =  ------ + ------ + ------
%         dxdx     dydy     dzdz
%                                                                                       2
%                  2  d   /  2  du  \                 d    /              du    \      d u
%     =  sin(theta)  ---- | r  ----  | + sin(theta) ------ | sin(theta) ------- | + ----------,
%                     dr  \     dr  /               dtheta \            dtheta  /    dphidphi
% where r = sqrt(x^2 + y^2 + z^2) is the radial distance from the origin,
% theta = acos(z/r) is the (geocentric) colatitude on (0, pi), and
% phi = atan2(y, x) is the (geocentric) longitude on (0, 2pi). These
% definitions of spherical coordinates are used with the physics convention
% (spherical coordinates defined with mathematics and physics conventions
% interchange the role of theta and phi). The series solution of the
% unknown function u in Laplace's equation is written
%                       inf    n
%                       ___   ___   /     n       -1-n \ /                               \                 
%                       \     \     | A  r  + B  r     | | C  cos(m phi) + S  sin(m phi) |  P  (cos(theta)),
%  u(r, theta, phi) =   /__   /__   \  n       n       / \  nm              nm           /   nm            
%                      n = 0  m = 0
% where An, Bn, Cnm, and Snm are the constants defined by initial and
% boundary conditions imposed on u and Pnm is the associated Legendre
% polynomial of degree n and order m.
% 
%    Inputs:
% 
%                 N - Provides a list of degrees for each individiual
%                     harmonic coefficient being evaluated.
%                     Size: L-by-1 (vector)*
%                     Units: - (N/A)
% 
%                 M - Provides a list of orders for each individiual
%                     harmonic coefficient being evaluated.
%                     Size: L-by-1 (vector)*
%                     Units: - (N/A)
% 
%                 C - Cosine-series harmonic coefficients to be evaluated.
%                     Varies with both degree n and order m.
%                     Size: L-by-1 (vector)*
%                     Units: ?
% 
%                 S - Sine-series harmonic coefficients to be evaluated.
%                     Varies with both degree n and order m.
%                     Size: L-by-1 (vector)*
%                     Units: ?
% 
%         longitude - Longitude (geocentric) corresponding with standard
%                     spherical coordinates, in which, given the
%                     coordinates of a position (x, y, z) in the Cartesian
%                     frame, the longitude is obtained according to
%                                  longitude = atan2(y, x).
%                     Thus, the longitude represents the angle made between
%                     the +x axis and the projection of the position into
%                     the x-y plane (precisely (x, y, 0)).
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%          latitude - Latitude (geocentric) corresponding with standard
%                     spherical coordinates, in which, given the
%                     coordinates of a position (x, y, z) in the Cartesian
%                     frame, the latitude is obtained according to
%                               latitude = asin(z / r).
%                     Note that the latitude is NOT the same as the
%                     colatitude (theta, represented above). They are
%                     related precisely by
%                              latitude + colatitude = pi/2.
%                     Thus, any quantity theta in the solution is replaced
%                     by the quantity (pi/2 - latitude). Further note that
%                            cos(colatitude) = sin(latitude)
%                     under this definition.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%            radius - Optional(!) Radial distance corresponding with
%                     standard spherical coordinates, in which, given the
%                     coordinates of a position (x, y, z) in the Cartesian
%                     frame, the radius is obtained according to
%                              radius = sqrt(x^2 + y^2 + z^2).
%                     Size: 1-by-1 (scalar)
%                     Units: ? (length)
% 
%                 A - Optional(!) Euler-series harmonic coefficients to be
%                     evaluated with the r^n term.
%                     Size: L-by-1 (vector)*
%                     Units: ?
% 
%                 B - Optional(!) Euler-series harmonic coefficients to be
%                     evaluated with the r^n term.
%                     Size: L-by-1 (vector)*
%                     Units: ?
% 
%    Outputs:
%               SHS - (S)pherical (H)armonic (S)eries evaluated according
%                     to the provided solution to Laplace's equation above
%                     using the harmonic coefficients An, Bn, Cnm, and Snm
%                     of up to degree Nmax and order Mmax, where Nmax and
%                     Mmax are respectively written max(n) and max(m).
%                     Size: 1-by-1 (scalar)
%                     Units: ?
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2,
%                     where (N, M) are the maximum degree and order of the
%                     model respectively.
% 

% Check that a valid combination of inputs has been provided
if (nargin ~= 6 && nargin ~= 9)
    error("Invalid amount of arguments.")
end
% Check size of coefficients provided to ensure that they are they same
checkCoefficients(C, S);
if (nargin == 9)
    checkCoefficients(A, B);
    checkCoefficients(A, C);
end

% Specify the unspecified arguments so that they disappear from the series

% Precalculate the argument of the associated Legendre polynomials
sinLatitude = sin(latitude);
% Begin recurrence relations trigonometric functions
cosLon = cos(longitude);
sinLon = sin(longitude);

% Get maximum degree and order of the series
[Nmax, Mmax] = deal(max(N), max(M));

% Track linear index within the summation
nmidx = 1;

% Begin summation
SHS = 0;
for nn = 0:Nmax
    % Precompute all necessary Legendre polynomials for this degree
    PnmAll = computefnLegendre(nn, sinLatitude);
    % Precalculate terms that only depend on the degree
    Anrn_plus_Bnrn1minusn = 1;
    if (nargin == 9)
        Anrn_plus_Bnrn1minusn = A(nmidx)*rn + B(nmidx)*(radius^-1)*(rn^-1);
    end
    sumOfOrders = 0;
    for mm = 0:min(nn, Mmax)
        % Increase efficiency by utilizing recurrence relations on
        % trigonometric functions
        if (mm > 0)
            cosmLon = cosmLonNext;
            sinmLon = sinmLonNext;
        else
            cosmLon = 1;
            sinmLon = 0;
        end
        % Obtain these coefficients of degree n and order m
        Cnm = C(nmidx);
        Snm = S(nmidx);
        Pnm = PnmAll(1+mm);
        sumOfOrders = sumOfOrders + (Cnm*cosmLon + Snm*sinmLon)*Pnm;
        % Recurrence relations
        cosmLonNext = cosLon*cosmLon - sinLon*sinmLon;
        sinmLonNext = sinLon*cosmLon + cosLon*sinmLon;
        % Increase count
        nmidx = nmidx + 1;
    end
    SHS = SHS + Anrn_plus_Bnrn1minusn*sumOfOrders;
end