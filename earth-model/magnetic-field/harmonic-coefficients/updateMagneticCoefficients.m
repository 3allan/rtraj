function [gnm, hnm] = updateMagneticCoefficients(t, gnmo, hnmo, dgnmdt, dhnmdt)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Update nominal gravitational coefficients to account for secular changes
% since 2000.0 and perturbations due to the effects of the solid Earth
% tides and ocean tides due to the ephemeris of the sun and moon.
% 
%    Inputs:
% 
%                 t - Time.
%                     Size: 1-by-1 (vector)*
%                     Units: JD UTC (Julian date in UTC)
% 
%               gnm - Cosine-series harmonic coefficients of the magnetic
%                     model of degree 12 and order 12 as defined at the epoch
%                     2020.0. These gnm are not normalized and carry
%                     dimensions of nT.
%                     Size: L-by-1 (vector)*
%                     Units: nT (nano-Tesla)
% 
%               hnm - Sine-series harmonic coefficients of the magnetic
%                     model of degree 12 and order 12 as defined at the epoch
%                     2020.0. These hnm are not normalized and carry
%                     dimensions of nT.
%                     Size: L-by-1 (vector)*
%                     Units: nT (nano-Tesla)
% 
%            dgnmdt - Constant secular rates of the cosine harmonic
%                     coefficients of the magnetic field model, gnm. They
%                     are also not normalized and carry units of nT/year.
%                     Size: L-by-1 (vector)*
%                     Units: nT/yr (nano-Tesla per year)
% 
%            dhnmdt - Constant secular rates of the cosine harmonic
%                     coefficients of the magnetic field model, gnm. They
%                     are also not normalized and carry units of nT/year.
%                     Size: L-by-1 (vector)*
%                     Units: nT/yr (nano-Tesla per year)
% 
%    Outputs:
% 
%               gnm - Corrected cosine-series magnetic harmonic coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: nT (nano-Tesla)
% 
%               hnm - Corrected sine-series magnetic harmonic coefficients.
%                     Size: L-by-1 (vector)*
%                     Units: nT (nano-Tesla)
% 
%                   * L = N - 2 - M (M - 2N - 1) / 2,    (N = M = 12).
% 

% Display
disp("Updating magnetic model...")

% Initialize output
gnm = gnmo;
hnm = hnmo;

% Obtain Julian date in UTC of 2020.0
JD2020UTC = juliandate([2020, 1, 1, 0, 0, 0]);

% Obtain Julian years since 2020.0
tSince2020JYUTC = (t - JD2020UTC) / 365.25;

% Apply corrections to the magnetic coefficients through simple linear
% variation in time since the epoch 2020.0
gnm = gnm + dgnmdt*tSince2020JYUTC;
hnm = hnm + dhnmdt*tSince2020JYUTC;

% Check nonzero elements again
checkCoefficients(gnm, hnm, "Magnetic", gnmo, hnmo)

% Display
fprintf("Updated magnetic model\n\n")
