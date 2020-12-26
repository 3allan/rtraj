function [g_sphr] = GradGravityPotential(GM, Req, ...
                                         n, m, C, S, ...
                                         r, colatitude, longitude)
% 
% Matt Werner (m.werner@vt.edu) - Dec 22, 2020
% 
% Calculate the gradient of the gravitational parameter in exact,
% pointwise fashion by analytically evaluating the gradient at the
% corresponding position indicated by the spherical coordinates r,
% colatitude (theta), and longitude (lambda) with respect to the ECF frame.
% The ECF frame provides a set of coordinates fixed to the rotating Earth,
% so the spherical coordinates here are essentially stationary with respect 
% to proper time. To be consistent with the coordinate system of the
% inputs, the result is returned in the same ECF spherical coordinate
% system. The result of this function is NOT the full gravitational 
% acceleration at the specified point since centrifugal acceleration due to
% the Earth's rotation must also be included.
% 
%    Inputs:
% 
% 
%                GM - Central-force gravitational parameter of the Earth
%                     consistent with the corresponding gravitational
%                     potential model defined by the harmonic coefficients
%                     C and S. The gravitational parameter is the
%                     standard product between the Earth's mass (M) and
%                     Newton's gravitational constant (G).
%                     Size: 1-by-1 (scalar)
%                     Units: m3/s2 (cubic meters per squared seconds)
% 
%               Req - Equatorial radius (semimajor axis) of the ellipsoid
%                     consistent with the corresponding gravitational
%                     potential model defined by the harmonic coefficients
%                     C and S.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 n - Degrees of the harmonic coefficients to be used in
%                     the spherical harmonic series expression of the
%                     gravitational potential.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                 m - Orders of the harmonic coefficients to be used in
%                     the spherical harmonic series expression of the
%                     gravitational potential.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                 C - Cosine-series harmonic coefficients to be used in the
%                     spherical harmonic series expression of the
%                     gravitational potential. The coefficients are defined
%                     up to a maximum degree N and order M corresponding
%                     with the largest values in the specified degrees n
%                     and orders m.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                 S - Sine-series harmonic coefficients to be used in the
%                     spherical harmonic series expression of the
%                     gravitational potential. The coefficients are defined
%                     up to a maximum degree N and order M corresponding
%                     with the largest values in the specified degrees n
%                     and orders m.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                 r - Radial distance of the current position from Earth's
%                     center (the origin of the ECF frame). In standard
%                     spherical coordinates, this quantity is given the
%                     symbol 'r'. For this model, appropriate values of r
%                     must officially abide by 
%                                       Req < r < inf, 
%                     but r slightly less than Req still leads to
%                     convergence of the series while incurring minimal
%                     errors.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%        colatitude - Colatitude of the current position relative to the
%                     ECF frame. This quantity is NOT the geocentric
%                     latitude, but rather the geocentric colatitude, which
%                     is the angle measured relative to deflection from the
%                     3-axis passing through the North pole. In standard
%                     spherical coordinates (used in physics), this
%                     quantity is given the symbol 'theta'. The colatitude
%                     must satisfy the inequality
%                                  0 <= colatitude <= pi,
%                     where colatitudes correspond to 
%                     1.    0 - a position above the North pole
%                     2. pi/2 - a position along the equatorial plane
%                     3.   pi - a position above the South pole.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%         longitude - Longitude of the current position relative to the ECF
%                     frame. In standard spherical coordinates (used in
%                     physics), this quantity is given the symbol 'phi' (or
%                     typically 'lambda' for applications directly related
%                     to physical longitudes related to the earth). The
%                     longitude must satisfy the inequality
%                                 -pi <= longitude < +pi,
%                     where longitudes correspond to
%                     1.     0 - a position along the same meridian arc as
%                               Greenwich, England
%                     2.  pi/2 - a position in the meridian arc 90 degrees
%                               east of Greenwich, England
%                     3. -pi/2 - a position in the meridian arc 90 degrees
%                               west of Greenwich, England.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    Outputs:
% 
%            g_sphr - The acceleration due to the gravitational attraction
%                     between the Earth and a small mass. This quantity is
%                     NOT the full gravitational acceleration since the
%                     effect of Earth's rotation is not taken into account
%                     here. It is given as a vector expressed in ECF
%                     spherical coordinates as (r, colatitude, longitude).
%                     Size: 3-by-1 (scalar)
%                     Units: m/s2 (meters per squared second)
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2,
%                     where (N, M) are the maximum degree and order of the
%                     model respectively.
% 

% Define spherical coordinate tuple and their appropriate ranges
sphrcoords = [r; colatitude; longitude];
sphrcoordsI = [0.95*Req, inf; -pi, pi; 0, 2*pi];

% Check spherical coordinate inputs
checkxInInterval(sphrcoords, sphrcoordsI)

% Get maximum degree and order
[N, M] = deal(max(n), max(m));

% Calculate the argument of the Legendre polynomials
coscolat = cos(colatitude);

% Precalculate some repeating terms to increase efficiency
GM_over_r2 = GM/r^2;
Req_over_r = Req/r;
Reqn_over_rn = Req_over_r;
Pnmcoscolat = computefnLegendre(2, coscolat);
cotcolat = cot(colatitude);
csccolat = csc(colatitude);

% Begin recurrence relations on trigonometric functions
cosLon = cos(longitude);
sinLon = sin(longitude);

% Start a counter to track the current position in C and S
nmidx = 1;

% Begin evaluating the gradient
g_sphr = zeros(3, 1);
for nn = 2:N
    % Reset inside summations for r, colat, and longitude directions
    orderSeries_r = 0;
    orderSeries_colat = 0;
    orderSeries_longitude = 0;
    % Increment the power of (Req/r) to the current value of nn
    Reqn_over_rn = Reqn_over_rn*Req_over_r;
    % Assign value for the next degree
    nnp1 = nn + 1;
    % Compute next degree of Legendre polynomials
    Pnp1mcoscolat = computefnLegendre(nnp1, coscolat);
    for mm = 0:min(nn, M)
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
        
        % Assign value for the next order
        mmp1 = mm + 1;
        
        % Assign this degree (and +1) and order's Legendre polynomial
        Pnm = Pnmcoscolat(mmp1);
        Pnp1m = Pnp1mcoscolat(mmp1);
        % Precompute the standard longitudinal variation of the series
        CnmCosmLon_plus_SnmSinmLon = Cnm*cosmLon + Snm*sinmLon;
        % Compute the derivative of Pnm(cos(x)) with respect to x using
        % this degree and the next degree evaluated at the same x and
        % order m with its leading negative sign disposed away (to be added
        % back in at the very end of this evaluation)
        negative_dPnmcoscolatdcolat = (nnp1 - mm)*csccolat*Pnp1m + nnp1*cotcolat*Pnm;
        
        % Compute the inner summations dependent upon the order m
        orderSeries_r = orderSeries_r + CnmCosmLon_plus_SnmSinmLon*Pnm;
        orderSeries_colat = orderSeries_colat + CnmCosmLon_plus_SnmSinmLon*negative_dPnmcoscolatdcolat;
        orderSeries_longitude = orderSeries_longitude + mm*(Cnm*sinmLon - Snm*cosmLon)*Pnm;
        
        % Recurrence relations
        cosmLonNext = cosLon*cosmLon - sinLon*sinmLon;
        sinmLonNext = sinLon*cosmLon + cosLon*sinmLon;
        % Increment counter
        nmidx = nmidx + 1;
    end
    % Apply degree-only multiplicative factors to the inside summations
    orderSeries_r = nnp1*Reqn_over_rn*orderSeries_r;
    orderSeries_colat = Reqn_over_rn*orderSeries_colat;
    orderSeries_longitude = Reqn_over_rn*orderSeries_longitude;
    
    % Add contribution of this degree series term to the total
    g_sphr = g_sphr + [orderSeries_r; orderSeries_colat; orderSeries_longitude];
    
    % Switch over the Legendre polynomials
    Pnmcoscolat = Pnp1mcoscolat;
end

% Account for Newtonian point-mass radial acceleration
g_sphr(1) = g_sphr(1) + 1;
% Apply cosecant term to the longitudinal acceleration
g_sphr(3) = g_sphr(3)*csccolat;
% Apply scaling and units to the gravitational acceleration due solely to
% massive attraction
g_sphr = -GM_over_r2*g_sphr;