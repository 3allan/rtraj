function JDy = convertJulianDayClock(JDx, xTimeSystem, yTimeSystem)
% 
% Matt Werner (m.werner@vt.edu) Dec. 1, 2020
% 
% Convert time systems used within the Julian day JDx. Relevant time
% systems are UT1, TAI, UTC, TT, and ET. These time systems are defined:
% 
% UT1 (Universal Time 1): Equivalent to Greenwich Mean Time (GMT) which 
% is the local time in Greenwich, England. UT1 operates on the constant
% 24-hour mean solar day, which does not account for the nonuniformity of
% Earth's rotation.
% 
% TAI (International Atomic Time): TAI is a realization for the statistical
% measure of time using about 200 atomic clocks tracking the amount of
% seconds past Jan 1, 1958 00:00:000 UT1. The official SI unit of 1 second 
% is precisely the amount of time it takes the Cesium-133 isotope at 
% absolute zero temperature to undergo exactly 9,192,631,770% periods of 
% radiation produced by transitions between the two hyperfine levels of 
% its unperturbed ground state, which serves as the basis for atomic
% clocks.
% 
% UTC (Coordinated Universal Time): Provides the link between UT1 (common
% time-keeping) with atomic time TAI (precise time-keeping) through the 
% introduction of leap seconds. A leap second is added whenever the time
% difference between UT1 and UTC exceeds about 0.6 seconds. A permanent,
% fixed amount of time separates UTC from TAI (exactly 10 seconds) along
% with the addition of leap seconds introduced into calendars.
% 
% TT (Terrestrial Time): Defined by NASA's Navigation and Ancillary 
% Information Facility (NAIF) to be the "proper" time on Earth at sea
% level. TT is offset from TAI by a permanent, fixed amount of exactly
% 32.184 seconds.
% 
% ET (Ephemeris Time): The mathematical idealization of the notion of time
% when writing equations of motion. It differs from TT by never more than
% 1.66 milliseconds and depends on the Earth's position (eccentric anomaly)
% about the sun according to Kepler's law and general relativity.
% 
%    Inputs:
% 
%               JDx - Julian date of a time t with respect to the time
%                     system x which is to be converted into the time
%                     system y. Here, x is usually any of UT1, TAI, UTC,
%                     TT, or ET.
%                     Size: 1-by-1 (scalar)
%                     Units: JD x (Julian date in the time system x)
% 
%   (x/y)TimeSystem - Time system of the provided Julian date (x) and
%                     desired time system to which JDx shall be converted
%                     (y) respectively.
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%                    Permissible options are:
%                         "UT1" - Common time-keeping time system that
%                                 operates on the constant 24-hour mean
%                                 solar day clock. UT1 is local to 
%                                 Greenwich, England.
% 
%                         "TAI" - Precise time-keeping time system that
%                                 operates according to a statistical
%                                 measure of about 200 atomic clocks
%                                 located in various labs throughout the
%                                 world.
% 
%                         "UTC" - Realization of connecting common time
%                                 (UT1) to atomic time (TAI) through the
%                                 introduction of leap seconds. TAI and UTC
%                                 are always kept within 1 second of each
%                                 other at most, though usually a leap
%                                 second is added when the deviation
%                                 between TAI and UTC reaches about 0.6
%                                 seconds. UTC and TAI are related simply
%                                   TAI = UTC + 10 + # of leap seconds,
%                                 where 10 is in seconds and the number of
%                                 leap seconds is the net count of leap
%                                 seconds added to UTC since the first was
%                                 added on Dec 31, 1971 23:59:60 UTC. This
%                                 formula relating UTC with TAI is valid
%                                 for all times after Jan 1, 1972 00:00:000
%                                 UTC.
% 
%                          "TT" - Recognized by NAIF to be the proper time
%                                 on Earth at sea level. TT and TAI are
%                                 both atomic time systems separated by a
%                                 permanent, fixed amount of 32.184
%                                 seconds with TT ahead of TAI. Thus, their
%                                 relation follows.
%                                            TT = TAI + 32.184,
%                                 where, again, the quantity 32.184 is
%                                 measured in seconds. Thus, TT and UTC
%                                 also have a very simple relation with
%                                 each other which just involves a constant
%                                 offset and a number of leap seconds.
% 
%                          "ET" - Recognized by NAIF to be the mathematically 
%                                 ideal notion of time for writing equations 
%                                 of motion. ET varies from TT by never
%                                 more than 1.66 milliseconds whose
%                                 variation is determined by the theory of
%                                 general relativity, though their relation
%                                 is still relatively straightforward.
%                                     ET = TT + 0.001657 sin(E + e sin E),
%                                 where 0.001657 is a quantity measured in
%                                 seconds, E is Earth's eccentric anomaly
%                                 about the sun, and e is Earth's orbital
%                                 eccentricity about the sun.
% 
%    Outputs:
% 
%               JDy - Julian date of a time t with respect to the time
%                     system y which was converted from the time system x.
%                     Here, y is usually any of UT1, TAI, UTC, TT, or ET.
%                     Size: 1-by-1 (scalar)
%                     Units: JD y (Julian date in the time system y)
% 

% Check that string inputs are strings (or characters)
checkInput(xTimeSystem, yTimeSystem);
% Convert to uppercase
[xTimeSystem, yTimeSystem] = deal(upper(xTimeSystem), upper(yTimeSystem));

% Convert time systems by adding or subtracting the correct amount of
% seconds from JDx.
switch xTimeSystem
    case {"UT1", "UTC"}
        switch yTimeSystem
            case {"UT1", "UTC"}
                JDy = JDx;
            case "TAI"
                % Get leap seconds
                leapseconds = getLeapSeconds;
                % Identify UT1 with UTC
                JDy = JDx + (10 + leapseconds)/86400;
            case "TT"
                % Copy UT1 to TAI and add it with 32.184 seconds to convert
                % to TT.
                % Get leap seconds
                leapseconds = getLeapSeconds;
                % Identify UT1 with UTC
                JDy = JDx + (42.184 + leapseconds)/86400;
            otherwise
                error("Please provide a valid time system to convert Julian dates into.")
        end
    case "TAI"
        switch yTimeSystem
            case {"UT1", "UTC"}
                % Get leap seconds
                leapseconds = getLeapSeconds;
                JDy = JDx - (10 + leapseconds)/86400;
            case "TAI"
                JDy = JDx;
            case "TT"
                JDy = JDx + 32.184/86400;
            otherwise
                error("Please provide a valid time system to convert Julian dates into.")
        end
    case "TT"
        switch yTimeSystem
            case {"UT1", "UTC"}
                % Get leap seconds
                leapseconds = getLeapSeconds;
                JDy = JDx - (42.184 + leapseconds)/86400;
            case "TAI"
                JDy = JDx - 32.184/86400;
            case "TT"
                JDy = JDx;
            otherwise
                error("Please provide a valid time system to convert Julian dates into.")
        end
    otherwise
        error("Please provide a valid time system to convert Julian dates from.")
end