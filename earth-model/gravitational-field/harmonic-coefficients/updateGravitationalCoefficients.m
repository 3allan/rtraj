function [Cnm, Snm] = updateGravitationalCoefficients(t, n, m, Cnmo, Snmo)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Update nominal gravitational coefficients to account for secular changes
% since 2000.0 and perturbations due to the effects of the solid Earth
% tides and ocean tides due to the ephemeris of the sun and moon.
% 
%    Inputs:
% 
%                 t - Current time
%                     Size: 1-by-1 (scalar)
%                     Units: JD UTC (Julian date in UTC)
% 
%                 n - List of degrees at which each harmonic coefficient
%                     Cnm and Snm are evaluated.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                 m - List of orders at which each harmonic coefficient
%                     Cnm and Snm are evaluated.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%              Cnmo - Original cosine-series harmonic coefficients of the
%                     gravity model of degree N and order M as defined at
%                     the epoch 2008.0 except C20, C30, and C40 which are 
%                     defined at 2000.0 (J2000). These Cnm are fully-
%                     normalized, non-dimensional, nominal (unperturbed)
%                     coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%              Snmo - Original sine-series harmonic coefficients of the
%                     gravity model of degree N and order M as defined at
%                     the epoch 2008.0. These Snm are fully-normalized, 
%                     non-dimensional, nominal (unperturbed) coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%    Outputs:
% 
%               Cnm - Corrected cosine-series fully-normalized, non-
%                     dimensional gravitational harmonic coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%               Snm - Corrected sine-series fully-normalized, non-
%                     dimensional gravitational harmonic coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2,
%                     where (N, M) are the maximum degree and order of the
%                     model respectively.
% 

%% Coefficients for correction
% Secular rates of fully-normalized nondimensional zonal harmonic
% coefficients C20, C30, and C40 as of J2000.
[dC20dt, dC30dt, dC40dt] = deal(11.6e-12, 4.9e-12, 4.7e-12);

% Nominal Love numbers for the solid Earth tide external potential
% (elastic)
ek(2,0+1) = 0.29525;
ek(2,1+1) = 0.29470;
ek(2,2+1) = 0.29801;
ek(3,0+1) = 0.093;
ek(3,1+1) = 0.093;
ek(3,2+1) = 0.093;
ek(3,3+1) = 0.094;
% plus variant (knm+)
ekp(2,0+1) = -0.00087;
ekp(2,1+1) = -0.00079;
ekp(2,2+1) = -0.00057;
% (anelastic)
ak(2,0+1) = 0.30190 + 0i;
ak(2,1+1) = 0.29830 - 0.00144i;
ak(2,2+1) = 0.30102 - 0.00130i;
% plus variant (knm+)
akp(2,0+1) = -0.00089;
akp(2,1+1) = -0.00080;
akp(2,2+1) = -0.00057;

WeightsAndipop{1} = ...
  [55.565,  0.00221,  0,  0,  0,  0,  1,  0.01347, -0.00541, 16.6, -6.7 
   55.575,  0.00441,  0,  0,  0,  0,  2,  0.01124, -0.00488, -0.1,  0.1 
   56.554,  0.04107,  0, -1,  0,  0,  0,  0.00547, -0.00349, -1.2,  0.8 
   57.555,  0.08214,  0,  0, -2,  2, -2,  0.00403, -0.00315, -5.5,  4.3 
   57.565,  0.08434,  0,  0, -2,  2, -1,  0.00398, -0.00313,  0.1, -0.1 
   58.554,  0.12320,  0, -1, -2,  2, -2,  0.00326, -0.00296, -0.3,  0.2 
   63.655,  0.47152,  1,  0,  0, -2,  0,  0.00101, -0.00242, -0.3,  0.7 
   65.445,  0.54217, -1,  0,  0,  0, -1,  0.00080, -0.00237,  0.1, -0.2 
   65.445,  0.54438, -1,  0,  0,  0,  0,  0.00080, -0.00237, -1.2,  3.7 
   65.465,  0.54658, -1,  0,  0,  0,  1,  0.00079, -0.00237,  0.1, -0.2 
   65.655,  0.55366,  1,  0, -2,  0, -2, -0.00077, -0.00236,  0.1, -0.2 
   73.555,  1.01590,  0,  0,  0, -2,  0, -0.00009, -0.00216,  0.0,  0.6 
   75.355,  1.08875, -2,  0,  0,  0,  0, -0.00018, -0.00213,  0.0,  0.3 
   75.555,  1.09804,  0,  0, -2,  0, -2, -0.00019, -0.00213,  0.6,  6.3 
   75.565,  1.10024,  0,  0, -2,  0, -1, -0.00019, -0.00213,  0.2,  2.6 
   75.575,  1.10245,  0,  0, -2,  0,  0, -0.00019, -0.00213,  0.0,  0.2 
   83.655,  1.56956,  1,  0, -2, -2, -2, -0.00065, -0.00202,  0.1,  0.2 
   85.455,  1.64241, -1,  0, -2,  0, -2, -0.00071, -0.00201,  0.4,  1.1 
   85.465,  1.64462, -1,  0, -2,  0, -1, -0.00071, -0.00201,  0.2,  0.5 
   93.555,  2.11394,  0,  0, -2, -2, -2, -0.00102, -0.00193,  0.1,  0.2 
   95.355,  2.18679, -2,  0, -2,  0, -2, -0.00106, -0.00192,  0.1,  0.1];

WeightsAndipop{2} = ...
 [125.755, 12.85429,  2,  0,  2,  0,  2,   -29,       3,     -0.1,  0.0 
  127.555, 12.92714,  0,  0,  2,  2,  2,   -30,       3,     -0.1,  0.0 
  135.645, 13.39645,  1,  0,  2,  0,  1,   -45,       5,     -0.1,  0.0 
  135.655, 13.39866,  1,  0,  2,  0,  2,   -46,       5,     -0.7,  0.1 
  137.455, 13.47151, -1,  0,  2,  2,  2,   -49,       5,     -0.1,  0.0 
  145.545, 13.94083,  0,  0,  2,  0,  1,   -82,       7,     -1.3,  0.1 
  145.555, 13.94303,  0,  0,  2,  0,  2,   -83,       7,     -6.8,  0.6 
  147.555, 14.02517,  0,  0,  0,  2,  0,   -91,       9,      0.1,  0.0 
  153.655, 14.41456,  1,  0,  2, -2,  2,   -168,     14,      0.1,  0.0 
  155.445, 14.48520, -1,  0,  2,  0,  1,   -193,     16,      0.1,  0.0 
  155.455, 14.48741, -1,  0,  2,  0,  2,   -194,     16,      0.4,  0.0 
  155.655, 14.49669,  1,  0,  0,  0,  0,   -197,     16,      1.3, -0.1 
  155.665, 14.49890,  1,  0,  0,  0,  1,   -198,     16,      0.3,  0.0 
  157.455, 14.56955, -1,  0,  0,  2,  0,   -231,     18,      0.3,  0.0 
  157.465, 14.56176, -1,  0,  0,  2,  1,   -233,     18,      0.1,  0.0 
  162.556, 14.91787,  0,  1,  2, -2,  2,   -834,     58,     -1.9,  0.1 
  163.545, 14.95673,  0,  0,  2, -2,  1,   -1117,    76,      0.5,  0.0 
  163.555, 14.95893,  0,  0,  2, -2,  2,   -1138,    77,    -43.4,  2.9 
  164.554, 15.00000,  0, -1,  2, -2,  2,   -1764,   104,      0.6,  0.0 
  164.556, 15.00000,  0,  1,  0,  0,  0,   -1764,   104,      1.6, -0.1 
  165.345, 15.02958, -2,  0,  2,  0,  1,   -3048,    92,      0.1,  0.0 
  165.535, 15.03665,  0,  0,  0,  0, -2,   -3630,   195,      0.1,  0.0 
  165.545, 15.03886,  0,  0,  0,  0, -1,   -3845,   229,     -8.8,  0.5 
  165.555, 15.07107,  0,  0,  0,  0,  0,   -4084,   262,    470.9, -30.2 
  165.565, 15.04328,  0,  0,  0,  0,  1,   -4355,   197,     68.1, -4.6 
  165.575, 15.04558,  0,  0,  0,  0,  2,   -4665,   334,     -1.6,  0.1 
  166.455, 15.07749, -1,  0,  0,  1,  0,   85693, 21013,      0.1,  0.0 
  166.544, 15.07993,  0, -1,  0,  0, -1,   35203,  2084,     -0.1,  0.0 
  166.554, 15.08214,  0, -1,  0,  0,  0,   22794,   358,    -20.6, -0.3 
  166.556, 15.08214,  0,  1, -2,  2, -2,   22780,   358,      0.3,  0.0 
  166.564, 15.08434,  0, -1,  0,  0,  1,   16842,   -52,     -0.3,  0.0 
  167.355, 15.11392, -2,  0,  0,  2,  0,   3755,   -189,     -0.2,  0.0 
  167.365, 15.11613, -2,  0,  0,  2,  1,   3552,   -182,     -0.1,  0.0 
  167.555, 15.12321,  0,  0, -2,  2, -2,   3025,   -160,     -5.0,  0.3 
  167.565, 15.12542,  0,  0, -2,  2, -1,   2892,   -154,      0.2,  0.0 
  168.554, 15.16427,  0, -1, -2,  2, -2,   1638,   -93,      -0.2,  0.0 
  173.655, 15.51259,  1,  0,  0, -2,  0,   370,   -20,       -0.5,  0.0 
  173.665, 15.51480,  1,  0,  0, -2,  1,   369,   -20,       -0.1,  0.0 
  175.445, 15.58323, -1,  0,  0,  0, -1,   325,   -17,        0.1,  0.0 
  175.455, 15.58545, -1,  0,  0,  0,  0,   324,   -17,       -2.1,  0.1 
  175.465, 15.58765, -1,  0,  0,  0,  1,   323,   -16,       -0.4,  0.0 
  183.555, 16.05697,  0,  0,  0, -2,  0,   194,   -8,        -0.2,  0.0 
  185.355, 16.12989, -2,  0,  0,  0,  0,   185,   -7,        -0.1,  0.0 
  185.555, 16.13911,  0,  0, -2,  0, -2,   184,   -7,        -0.6,  0.0 
  185.565, 16.14131,  0,  0, -2,  0, -1,   184,   -7,        -0.4,  0.0 
  185.575, 16.14352,  0,  0, -2,  0,  0,   184,   -7,        -0.1,  0.0 
  195.455, 16.68348, -1,  0, -2,  0, -2,   141,   -4,        -0.1,  0.0 
  195.465, 16.68569, -1,  0, -2,  0, -1,   141,   -4,        -0.1,  0.0];

WeightsAndipop{3} = ...
 [245.655, 28.43973,  1,  0,  2,  0,  2,   0.00006, 0,       -0.3,  0.0 
  255.555, 28.98410,  0,  0,  2,  0,  2,   0.00004, 0,       -1.2,  0.0];

% Display
disp("Updating gravity model...")

% Initialize output
Cnm = Cnmo;
Snm = Snmo;

%% Correct for secular rate from J2000
% Get J2000 epoch
J2000 = getJulianDateofJ2000; % [Julian days]

% Calculate (Julian) days since J2000 (Jan 1, 2000 12:00:000 TT)
tSinceJ2000JD = t - J2000;
% Convert to (Julian) years since J2000
tSinceJ2000JY = tSinceJ2000JD/365.25;

% Pull out the nominal zonal coefficients (Cn0)
Cn0nom = Cnmo(m == 0);

% Correct the first 3 of them by applying a linear variation in their
% secular rates
numelCn0 = numel(Cn0nom);
dCn0dt = [dC20dt; dC30dt; dC40dt; zeros(numelCn0 - 3, 1)];
Cnm(m == 0) = Cn0nom + dCn0dt*tSinceJ2000JY;

%% Correct for the effect of solid Earth tides
% NASA JPL ephemerides of celestial bodies - module 432t released in April
% 2014. The ephemeris gives position (and velocity) of bodies relative to
% the International Celestial Reference Frame (v2.0 for this module)
% (ICRF). SPICE routines do not differentiate between ICRF and J2000 (ECI
% frame at epoch J2000) frames at all. Therefore, positions and angles
% gives directly from the ephemeris are given with respect to J2000 ECI
% frame. The Greenwich hour angle is therefore needed to transform these
% positions into the standard Earth-centered-fixed (ECF) frame.

% Get parameters of the earth
[mue, ~, ~, ~, ~, ~] = defineEllipsoidParameters("WGS84");
Re = 6379.194623; % Set radius of Earth for EGM2008 in km
% Get gravitational parameters of the sun and moon
[GMsun, GMmoon] = getGMsOfSunAndMoon;
% Get positions of sun and moon with respect to ICRF (~J2000 frame)
posLunaECI = planetEphemeris(t, 'Earth', 'Moon', '432t')';
posSolECI = planetEphemeris(t, 'Earth', 'Sun', '432t')';
% Get the distance (rotations leave distance invariant)
distLuna = norm(posLunaECI);
distSol = norm(posSolECI);
% Obtain rotation matrix from ECI to ECF coordinates at time t
Tecf_eci = getTransformationECI2ECFCoordinates(t);
% Transform positions of sun and moon to ECF frame
posLunaECF = Tecf_eci*posLunaECI;
posSolECF = Tecf_eci*posSolECI;
% Obtain the geocentric longitude and latitude of the moon
lonLunaECF = atan2(posLunaECF(2), posLunaECF(1));
latLunaECF = atan2(posLunaECF(3), sqrt(posLunaECF(1)^2 + posLunaECF(2)^2));
% Obtain the geocentric longitude and latitude of the sun
lonSolECF = atan2(posSolECF(2), posSolECF(1));
latSolECF = atan2(posSolECF(3), sqrt(posSolECF(1)^2 + posSolECF(2)^2));


% Begin calculating corrections
mu(2:3) = [GMmoon, GMsun]';
r(2:3, 1) = [distLuna, distSol]';
Phi(2:3, 1) = [latLunaECF, latSolECF]';
lambda(2:3, 1) = [lonLunaECF, lonSolECF]';
for nn = 2:3
    for mm = 0:nn
        sumEph = 0;
        for jj = 2:3
            % Place variables in ideal form to be indexed by j
            [muj, rj, Phij, lambdaj] = deal(mu(jj), r(jj), Phi(jj), lambda(jj));
            % Calculate the fully-normalized associated Legendre functions
            % for all 0 <= m <= n.
            PnmSbarsinPhij = sqrt(2*nn + 1) * legendre(nn, sin(Phij), 'sch');
            % Pick the one corresponding to m = m (1st element is m = 0)
            PnmbarsinPhij = PnmSbarsinPhij(mm + 1);
            % Form sum
            sumEph = sumEph + (muj/mue)*(Re/rj)^(nn + 1) * PnmbarsinPhij * exp(-1i*mm*lambdaj);
        end
        % Access the mm+1 element of knm since m = 0 isn't a valid
        % container in Matlab (they were all stored at index m + 1)
        knm = ek(nn, mm+1);
        if (nn == 2), knm = knm + ak(nn, mm+1); end
        DCnm_minus_iDSnm = (knm / 2*nn + 1)*sumEph;
        % Remove correction for the cosine-harmonic (2,0) term (set the
        % real part to nothing)
        if (nn == 2 && mm == 0)
            DCnm_minus_iDSnm = 0 + imag(DCnm_minus_iDSnm); 
        end
        % Apply corrections
        Cnm(n == nn & m == mm) = Cnm(n == nn & m == mm) + real(DCnm_minus_iDSnm);
        Snm(n == nn & m == mm) = Snm(n == nn & m == mm) - imag(DCnm_minus_iDSnm);
    end
end

% Apply a correction to the 4th degree coefficients as well induced by
% degree 2 tides
for mm = 0:2
    sumEph = 0;
    for jj = 2:3
        % Place variables in ideal form to be indexed by j
        [muj, rj, Phij, lambdaj] = deal(mu(jj), r(jj), Phi(jj), lambda(jj));
        % Calculate the fully-normalized associated Legendre functions
        % for all 0 <= m <= n.
        PnmSbarsinPhij = sqrt(5) * legendre(2, sin(Phij), 'sch');
        % Pick the one corresponding to m = m (1st element is m = 0)
        PnmbarsinPhij = PnmSbarsinPhij(mm + 1);
        % Form sum
        sumEph = sumEph + (muj/mue)*(Re/rj)^3 * PnmbarsinPhij * exp(-1i*mm*lambdaj);
    end
    % Sum up Love numbers of elastic and anelastic Earth
    knmp = ekp(2, mm+1) + akp(2, mm+1);
    DC4m_minus_iDS4m = (knmp / 5)*sumEph;
    % Apply corrections
    Cnm(n == 4 & m == mm) = Cnm(n == 4 & m == mm) + real(DC4m_minus_iDS4m);
    Snm(n == 4 & m == mm) = Snm(n == 4 & m == mm) - imag(DC4m_minus_iDS4m);
end
    
% ==== Something very wrong in this section === %
% Known: The multiplier (-1i)^mm is wrong, but there's more wrong with it
% % Apply step 2 correction
% 
% % Convert t to JC (TT) since J2000
% JDJ2000TT = getJulianDateofJ2000("TT");
% tJDTT = convertJulianDayClock(t, "UTC", "TT");
% tSinceJ2000JDTT = tJDTT - JDJ2000TT;
% tSinceJ2000JCTT = tSinceJ2000JDTT / 36525; % JC since J2000 (TT)
% % Assemble coefficients together for algebraic multiplication
% tvec = [1, tSinceJ2000JCTT, tSinceJ2000JCTT^2, tSinceJ2000JCTT^3, tSinceJ2000JCTT^4]';
% lvec = deg2rad([134.96340251, [1717915923.2178, 31.8792, 0.051635, -0.00024470]/3600]);
% lpvec = deg2rad([357.52910918, [129596581.0481, -0.5532, 0.000136, -0.00001149]/3600]);
% Fvec = deg2rad([93.27209062, [1739527262.8478, -12.7512, -0.001037, 0.00000417]/3600]);
% Dvec = deg2rad([297.85019547, [1602961601.2090, -6.3706, 0.006593, -0.00003169]/3600]);
% Ovec = deg2rad([125.04455501, [6962890.5431, 7.4722, 0.007702, -0.00005939]/3600]);
% % Compute the 5 Delaunay variables at time t in Dynamical time
% F = [lvec; lpvec; Fvec; Dvec; Ovec]*tvec;
% % Obtain the Greenwich Mean Sidereal Time (GMST)
% [~, GMST, ~] = getRotAngsfromJDUT1(t, nan);
% % Allocate space
% thetaf{1} = zeros(size(WeightsAndipop{1}, 1), 1);
% thetaf{2} = zeros(size(WeightsAndipop{2}, 1), 1);
% thetaf{3} = zeros(size(WeightsAndipop{3}, 1), 1);
% % Begin looping to make thetaf 
% for mm = 0:2
%     mmindx = mm + 1;
%     ip = WeightsAndipop{mmindx}(:, 10);
%     op = WeightsAndipop{mmindx}(:, 11);
%     Amp = ip + 1i*op;
%     for ii = 1:numel(thetaf{mmindx})
%         % Form row vector of weightings for this tidal frequency f
%         N = WeightsAndipop{mmindx}(ii, 3:7);
%         % Compute the angle theta for this tidal frequency f
%         thetaf{mmindx}(ii) = mm*(GMST + pi) - N*F;
%     end
%     % Start the step 2 correction now that we have everything for this m
%     DCnm_minus_iDSnm = (-1i)^mm * sum(Amp.*exp(1i*thetaf{mmindx}));
%     Cnm(n == 2 & m == mm) = Cnm(n == 2 & m == mm) + real(DCnm_minus_iDSnm);
%     Snm(n == 2 & m == mm) = Snm(n == 2 & m == mm) - imag(DCnm_minus_iDSnm);
% end
% ============================================= %

% Step 3 only applies if using the zero-tide gravitational model
% --> use tide-free for now

%% Correct for effects that ocean tides have on the gravitational field
% load("mats/EGM/coefficientmats/fes2004CnmSnm.mat", "fes2004CnmSnm", '-mat');
%  Match tidal frequencies used in thetaf with tidal frequencies in this .mat.

% Check how many coefficients have changed
checkCoefficients(Cnm, Snm, "Gravity", Cnmo, Snmo)

fprintf("Updated gravity model\n\n")
