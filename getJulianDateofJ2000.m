function J2000 = getJulianDateofJ2000(yTimeSystem)
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Get the Julian date (JD) of the epoch J2000 (Jan 1, 2000 12:00:000 GMT).
% Note that the Julian day is defined to have vanishing decimal part at
% noon in England - not midnight in England.
% 
%    Inputs:
%       yTimeSystem - Optional(!) input to specify which time system to
%                     convert the Julian date of J2000 into. The default
%                     time system is Terrestrial Time (TT) which is
%                     returned if no input is specified. The options are
%                     UT1, TAI, UTC, TT, and ET.
%                     Size: 1-by-1 (scalar)
%                     Units: JD UTC (Julian date in UTC)
% 
%    Outputs:
% 
%             J2000 - Julian date (JD) of the epoch J2000 in Terrestrial 
%                     time (TT). J2000 is defined to have occured precisely
%                     at the moment of Jan 1, 2000 at 12:00:000 TT - that 
%                     is, New Year's day at noon in Greenwich, England in 
%                     the year 2000.
%                     Size: 1-by-1 (scalar)
%                     Units: JD TT (default) (Julian date in TT or other)
% 

% Obtain Julian date of J2000 in TT by default
J2000 = 2451545.0; % TT
    
if (nargin == 1)
    J2000 = convertJulianDayClock(J2000, "TT", yTimeSystem);
end