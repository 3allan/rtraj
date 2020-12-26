function [n, m, Jn, Cnm, Snm, dCnm, dSnm] = loadGravitationalCoefficients(N, M, GMtype)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Load the gravitational coefficients defined by the Earth Gravity Model of
% 2008 (EGM2008) up to degree N and order M. The coefficients are loaded by
% default as the fully-normalized flavor and are completely unitless. The 
% epoch of the coefficients is precisely 2008.0 (Jan 1 2008). The
% uncertainties of each coefficient is also provided should they be desired.
% EGM2008 provides two (2) sets of gravitational coefficients - one
% pertaining to the tide-free Earth and the other pertaining to the
% zero-tide Earth
% 
% The tide-free gravitational model removes all direct and indirect tidal
% effects of the sun and moon.
% 
% The zero-tide gravitational model removes all permanent direct effects of
% the sun and moon, but retains the indirect effect related to the earth's 
% elastic deformation.
% 
% The only difference between the tide-free and zero-tide models lies in
% the expression of C20; the rest of the coefficients are identically the
% same.
% 
% If only N is specified, then the model is loaded as best as it can be to
% degree and order N. If M is zero, then all Snm are zero and the
% gravitational potential is independent of geodetic longitude, i.e. 
% equivalent to the normal gravity field. In this case, the Cnm and Jn are
% related.
% 
% EGM2008 is given in spherical coordinates relative to the standard Earth-
% Centered-Rotating (ECF) frame.
% 
%    Inputs:
% 
%                 N - Maximum desired degree of the gravity model. EGM2008
%                     provides a maximum degree of 2190. The gravity model
%                     is called "complete" such that it has a "full" set of
%                     coefficients defined up to degree 2159. "Full" here
%                     means that there are no missing orders m for a given
%                     degree n (a missing order m for a given degree n
%                     indicates that data are missing from the set).
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%                 M - Maximum desired order of the gravity model. EGM2008
%                     provides a maximum order of 2159. Orders vary for
%                     each degree n as 0 <= m <= n. Thus, EGM2008 is
%                     complete up to degree and order 2159. Additional
%                     coefficients are obtainable for degrees exceeding
%                     2159, but the coefficients will be zero (i.e. non-
%                     existent) for orders exceeding 2159.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%            GMtype - Name of model for the gravitational field of Earth.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                     "tide-free" - Provides the fully-normalized tide-free
%                                   unitless EGM2008 harmonic coefficients
%                                   Cnm and Snm at the epoch 2008.0 (Jan 1,
%                                   2008) except for the zonal coefficients 
%                                   C20, C30, and C40, which have epoch
%                                   2000.0 (J2000). These 3 coefficients
%                                   may be brought up to date by applying a
%                                   linear correction in time, where the
%                                   secular variation of each coefficient
%                                   is a known constant with units (1/year).
%
%                     "zero-tide" - Provides the fully-normalized zero-tide
%                                   unitless EGM2008 harmonic coefficients
%                                   Cnm and Snm at the epoch 2008.0 (Jan 1,
%                                   2008) except for the zonal coefficients 
%                                   C20, C30, and C40, which have epoch
%                                   2000.0 (J2000). These 3 coefficients
%                                   may be brought up to date by applying a
%                                   linear correction in time, where the
%                                   secular variation of each coefficient
%                                   is a known constant with units (1/year).
% 
%    Outputs:
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
%                Jn - Cosine-series harmonic coefficients of the gravity
%                     model of degree N and order 0 as defined at the epoch
%                     2008.0 except J2, J3, and J4 which are defined at
%                     2000.0 (J2000). The Jn series is closely related to
%                     the fully-normalized coefficients Cn0 by
%                                   Jn = -sqrt(2n + 1) Cn0,
%                     where it's again emphasized that Cn0 are the fully-
%                     normalized coefficients here. Thus, the prefactor is
%                     effectively un-normalizing the Cn0 coefficients.
%                     Size: (N-1)-by-1 (vector)
%                     Units: - (unitless)
% 
%               Cnm - Cosine-series harmonic coefficients of the gravity
%                     model of degree N and order M as defined at the epoch
%                     2008.0 except C20, C30, and C40 which are defined at
%                     2000.0 (J2000). These Cnm are fully-normalized, non-
%                     dimensional, nominal (unperturbed) coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%               Snm - Sine-series harmonic coefficients of the gravity
%                     model of degree N and order M as defined at the epoch
%                     2008.0 except C20, C30, and C40 which are defined at
%                     2000.0 (J2000). These Snm are fully-normalized, non-
%                     dimensional, nominal (unperturbed) coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%              dCnm - Uncertainties of the cosine-series harmonic 
%                     coefficients of the gravity model of degree N and 
%                     order M.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%              dSnm - Uncertainties of the sine-series harmonic 
%                     coefficients of the gravity model of degree N and 
%                     order M.
%                     Size: L-by-1 (vector)*
%                     Units: - (unitless)
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2
% 

% Display
disp("Loading gravity model...")

% Check if the input is of type string (or character))
checkInput(GMtype)
% Check if degree is specified too high
if (N > 2190)
    error("Cannot supply degrees N exceeding 2190.")
end
% Check if order is specified too high
if (M > N)
    error("Cannot supply orders M exceeding degree N.")
end
% Check if order is specified too high (again)
if (M > 2159)
    warning("Changing order M to maximum order (2159).")
    M = 2159;
end

% Get the coefficients
switch GMtype
    case "tide-free"
        load('mats/EGM/coefficientmats/EGM2008GravityCoeff_TideFree.mat', '-mat', 'EGM2008GravityCoeff_TideFree');
        G = EGM2008GravityCoeff_TideFree;
        clear EGM2008GravityCoeff_TideFree % Clear from memory (~118MB)
    case "zero-tide"
        load('mats/EGM/coefficientmats/EGM2008GravityCoeff_ZeroTide.mat', '-mat', 'EGM2008GravityCoeff_ZeroTide');
        G = EGM2008GravityCoeff_ZeroTide;
        clear EGM2008GravityCoeff_ZeroTide % Clear from memory (~118MB)
    otherwise
        error("Please indicate a valid variation of the gravity model to load.")
end

% Retrieve only the coefficients according to (N, M)
Gmaxrow = -2 + nchoosek(N + 1, 2) + M;
GNM = G(1:Gmaxrow, :); % Gravity potential coefficients

% Provide "corrections" by removing orders m > M
GNM(GNM(:, 2) > M, :) = [];

% Distribute results
[n, m, Cnm, Snm, dCnm, dSnm] = deal(GNM(:, 1), GNM(:, 2), GNM(:, 3), GNM(:, 4), GNM(:, 5), GNM(:, 6));

% Calculate Jn
nCn0 = GNM(GNM(:, 2) == 0, [1, 3]);
clear GNM
Jn = -sqrt(1 + 2*nCn0(:, 1)).*nCn0(:, 2);

% Display
L = N - 2 - M*(M - 2*N - 1)/2;
checkCoefficients(Cnm, Snm, "Gravity")
if (L == Gmaxrow || (N == 2190 && M == 2159))
    disp("Gravity model is full degree and order (degree " + N + " and order " +  M + ")")
end

% Display
fprintf("Loaded gravity model\n\n")