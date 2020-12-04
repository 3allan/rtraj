function [ERA, GMST, LST] = getRotAngsfromJDUT1(t, longitude)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Obtain expressions for the Earth Rotation Angle (ERA) and Greenwich Mean
% Sidereal Time (GMST) through the time of day expressed in Julian days in
% the UT1 (~UTC) time system.
% 
%    Inputs:
%                 t - Time
%                     Size: 1-by-1 (scalar)
%                     Units: JD UT1 (Julian date according to UT1)
% 
%         longitude - Optional(!) Geodetic longitude (or geocentric 
%                     longitude, they are the same in this case) of a 
%                     local location on Earth. The longitude is measured 
%                     relative to Greenwich, England, where a positive 
%                     longitude corresponds to the eastward direction and 
%                     a negative longitude corresponds to the westward 
%                     direction. Per standard definition of the longitude, 
%                     the domain varies between -180 < longitude < 180, 
%                     where longitude here is measured in degrees.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    Outputs:
%               ERA - Earth Rotation Angle (ERA) at the time t. The ERA is
%                     defined as the angle measured along the equator of
%                     the Celestial Intermediate Pole (CIP) between the
%                     Celestial Intermediate Origin (CIO) and Terrestrial
%                     Intermediate Origin (TIO). It follows a linear
%                     relation with UT1 as
%                               ERA = ERA0 + (dERA/dt)(t - t0),
%                     where ERA0 is the ERA at t = t0, dERA/dt is the rate
%                     of advance of ERA in rev/(UT1 day), and t0 = J2000 in
%                     JD UT1.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%              GMST - Greenwich Mean Sidereal Time (GMST) at the time t.
%                     The GMST is defined as the angle between the J2000
%                     vernal equinox (+x axis) and Greenwich, England. A
%                     positive angle corresponds with how far (in the
%                     angular sense) Greenwich, England has rotated away
%                     from the vernal equinox in the day. A modern
%                     expression for the GMST is given as
%                                   GMST = ERA + p(t),
%                     where p(t) is a quintic polynomial in time t as
%                     expressed in Julian days with respect to terrestrial
%                     time (TT).
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%               LST - Optional(!) Local Sidereal Time (LST) at the time t.
%                     The LMST is defined as the angle between the J2000
%                     vernal equinox (+x axis) to Greenwich, England, and
%                     then to a local location defined by a geodetic
%                     longitude measured positive eastward from Greenwich,
%                     England.The expression for the LMST is given as
%                             LMST = GMST + geodeticLongitude,
%                     where GMST is the Greenwich Mean Sidereal Time and
%                     geodeticLongitude is the geodetic (or geocentric)
%                     longitude of the local location.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 

% Prepare to calculate ERA
J2000TT = getJulianDateofJ2000;
J2000UT1 = getJulianDateofJ2000("UT1");
Tu = t - J2000UT1;
JDfracUT1 = mod(t, 1);

% Calculate Earth Rotation Angle (ERA) from IERS TN36 [rad]
ERA = mod(2*pi*(JDfracUT1 + 0.7790572732640 + 0.00273781191135448*Tu), 2*pi);

% Prepare to calculate GMST according to rtraj documentation (lost source?)
tTTJC = convertJulianDayClock(t - J2000UT1, "UT1", "TT") / 36525;
tTTJCvec = [1, tTTJC, tTTJC^2, tTTJC^3, tTTJC^4, tTTJC^5]';
GMSTtcoeffs = deg2rad([0.014506, 4612.156534, 1.3915817, -0.00000044, -0.000029956, -0.0000000368]/3600);
% Calculate Greenwich Mean Sidereal Time (GMST) [rad]
GMSTdoc = mod(ERA + GMSTtcoeffs*tTTJCvec, 2*pi);

% Calculate GMST from Vallado
TUT1 = (t - J2000TT) / 36525; % Convert to Julian centuries since J2000 (TT)
GMSTseconds = 67310.54841 + (876600*3600 + 8640184.812866)*TUT1 + 0.093104*TUT1^2 - 6.2e-6*TUT1^3;
GMSTseconds = mod(GMSTseconds, 86400); % Restrict domain back to 1 day (86400 seconds)
GMSTVallado = (GMSTseconds/240)*(pi/180); % Convert to radians: 1 second ~ 1/240 degrees

% Average the results since they're similar but not quite the same
GMST = (GMSTdoc + GMSTVallado) / 2;

% Check if a longitude was given as a second argument
if (nargin == 2)
    % Calculate Local Sidereal Team (LST) according to the given longitude
    LST = mod(GMST + longitude, 2*pi);
end