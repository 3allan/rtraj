% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Load the propulsion data, which requires specification of AT LEAST the
% thrust profile as a function of burn time. Additional specification of
% the mass flow rate is not required, but will be assumed to be a linear
% function of burn time varying between the stage's wet and dry masses if
% left unspecified. The specification of the chamber pressure profile is
% also not required. Providing the chamber pressure profile, however,
% yields extra information about isentropic flow properties of the exhaust
% through the nozzle.
% 

% SCRIPT - Creates flags to determine which models may be used again from
% the rtraj.m workspace (if already populated by a previous simulation) and
% which require reloading
refPropulsionCache

% Only load the propulsion profiles that are either not currently stored in
% the rtraj.m workspace or have changed from the previous simulation
% 
% Load thrust profile(s) from BurnSim csv file(s)
if (flags.load.ThrustProfile)
    for ii = loadTheseFTProfiles
        FTProfile{ii, 1} = readThrustProfile(FTprofile(ii));
    end
else
    disp("Hitting thrust profile cache...")
    fprintf("Loaded thrust profile\n\n")
end
% Obtain mass flow rate profile from BurnSim csv files
if (flags.load.MassFlowRateProfile)
    MFProfile = {stripcsv('csvs/S1SL.csv', [1, 6], 1); stripcsv('csvs/S235.csv', [1, 6], 1)};
    % Units of mass flow rate
    MFUnits = {["s", "lb/s"]; ["s", "lb/s"]};
else
    disp("Hitting mass flow rate cache...")
    fprintf("Loaded mass flow rate profile\n\n")
end
% Obtain chamber pressure profile from BurnSim csv files
if (flags.load.ChamberPressureProfile)
    PCProfile = {NaN; NaN};
    % Units of chamber pressure profile
    PCUnits = {["NaN", "NaN"]; ["NaN", "NaN"]};
else
    disp("Hitting chamber pressure cache...")
    fprintf("Loaded chamber pressure cache\n\n")
end