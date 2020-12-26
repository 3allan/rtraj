% The (R)ocket (Traj)ectory Program (rtraj)
% Matt Werner (m.werner@vt.edu)
clearvars

% Notes
% 1. Use a mono font to correctly display formatting of comments/tables
% 2. Disable wrap-lines in Command Window preferences to view tables easily
% 3. Automatically resets any settings changed by producing a run
% 4. Requires Aerospace toolbox

% Options
% The vector i:j with 0 < i < j < 9 means show phases i to j where
% Stage 1
% 1 - Rail clearance        3 - Separation Delay        5 - Parachute
% 2 - Burnout 1             4 - Apogee                  6 - Touchdown
% Stage 2
% 1 - Rail clearance        3 - Separation Delay        5 - Burnout 2        7 - Parachute
% 2 - Burnout 1             4 - Ignition Delay          6 - Apogee           8 - Touchdown
tablePhse = {1:6}; % i:j - Sequences to display in tables (associated with DispTable)
plotsPhse = {1:6}; % i:j - Sequences to display in plots (associated with Plots)
odestats_ = 'on'; % 'on'/'off' - Give ODE solver statistics
oderefine = 3; % Positive integer - Number of points to evaluate between time steps
odereltol = 1e-5; % Positive scalar - Relative error tolerance bound
odeabstol = 1e-7; % Positive scalar - Absolute error tolerance bound
odeinstep = 1e-3; % Positive scalar - Initial step size in time
odemxstep = inf; % Positive scalar - Maximum permissible time step (including inf)
HavePrtrb = false; % true/false - Consider perturbations to nominal flight
CloneData = false; % true/false - Save time by performing no additional calculations
ShowPlots = false; % true/false - Show plots
Runtime__ = true; % true/false - Show runtime after finishing
ODEtime__ = true; % true/false - Show solver times at each time step

%% Inputs
% Obtain actual time right now
tLocalNow = datetime('now', 'timezone', 'America/New_York');
% Obtain actual launch time as Julian date in UT1
[JDLaunch, tLaunchUTC, tLaunch] = defineLaunchTime('15-July-2091 13:53:45.234', 'd-MMMM-yyyy HH:mm:ss.SSS', '-05:00');
% Provide launch site location
LonLaunch = -119.0560102; % Longitude of launch site [deg]
LatLaunch = + 40.910733; % Geodetic latitude of launch site [deg]
% Provide launch site temperature and density
localTemp = 70; % Local surface temperature [F]
% Define launch tower length
towerSpan = convUnits(30, "feet", "meters"); % Length of launch tower [m]
% Define rocket parameters (lengths/distances)
numStages = 2; % Number of vehicle stages
LenRocket = [6.35; 3.71]; % Overall length of the rocket during each stage [m]
LenNoseCM = [4.15; 1.50]; % Distance from the nose cone's tip to the center of mass of the fully fueled rocket for each stage [m]
diaOuter_ = convUnits([8.50; 6.00], "inches", "meters"); % Rocket casing's outermost diameter for each stage [m]
diaThroat = convUnits([2.50; 1.60], "inches", "meters"); % Nozzle throat diameter [m]
diaExit__ = convUnits([6.01; 4.47], "inches", "meters"); % Nozzle exit diameter [m]
diaFlatDM = convUnits([36, 108; 36, 108], "inches", "meters"); % Flattened drogue & main parachute diameters (0 indicates no chute) [m]
% Define rocket parameters (propulsion)
massInit_ = [194.43; 65.33]; % Initial mass of the rocket at the beginning of each burn time [kg]
massMotor = [63.49; 23.58]; % Motor mass of each stage [kg]
burnTimes = [6.95; 5.09]; % Burn time of each stage (takes precedence over max time in FTProfile) [s]
FTProfile = {stripcsv('csvs/S1SL.csv', [1, 4], 1); stripcsv('csvs/S235.csv', [1, 4], 1)}; % Obtain thrust profile from BurnSim csv files
FTUnits__ = {["s", "lbf"]; ["s", "lbf"]}; % Units of thrust profile
MFProfile = {stripcsv('csvs/S1SL.csv', [1, 6], 1); stripcsv('csvs/S235.csv', [1, 6], 1)}; % Obtain mass flow rate profile from BurnSim csv files
MFUnits__ = {["s", "lb/s"]; ["s", "lb/s"]}; % Units of mass flow rate
PCProfile = {NaN; NaN}; % Obtain chamber pressure profile from BurnSim csv files
PCUnits__ = {["NaN", "NaN"]; ["NaN", "NaN"]}; % Units of chamber pressure profile
Yexhaust_ = [1.25; 1.25]; % Perfect ratio of specific heats for motor exhaust gases []
% Define rocket parameters (aerodynamics/geometry)
RasAeroCD = {csvread('csvs/H7512.csv', 1, 0); csvread('csvs/H7522.csv', 1, 0)}; % Obtain drag coeff. profile from RasAero csv files
% ... geometry here to define our own CD profile
% Define rocket parameters (staging behavior)
delaySep_ = [2; NaN]; % Time it takes for stages to separate after burn-out [s]
delayIgn_ = [0; NaN]; % Time it takes for stages to ignite after separation [s]
% Define rocket parameters (controlled descent)
deployAlt = convUnits([2000; 2000], "feet", "meters"); % Altitude of main parachute deployment [m]
% Constants and measured quantities relevant to ODE initial conditions
Rail2Rckt = diaOuter_(1)/2; % Perpendicular distance from rail to rocket's centerline [m]
Rail2Vert = 5; % Angle made between railing and vertical [deg]
East2DwnR = 30; % Angle made between launch (downrange) direction and due East [deg]
veps = 0; % Very small initial velocity [m/s]

% Obtain parameters that define an Earth model
[GM, Req, Rpo, f, e, w] = defineEllipsoidParameters('WGS84');
[nG, mG, ~, Cnm, Snm, ~, ~] = loadGravitationalCoefficients(1000, 1000, 'tide-free');
[Cnm, Snm] = updateGravitationalCoefficients(JDLaunch, nG, mG, Cnm, Snm);
[nM, mM, gnm, hnm, dgnmdt, dhnmdt] = loadMagneticCoefficients(12, 12);
[gnm, hnm] = updateMagneticCoefficients(JDLaunch, gnm, hnm, dgnmdt, dhnmdt); clear dgnmdt dhnmdt
[longitudes, geodeticLatitudes, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain('degrees');
atmosModel = "NRLMSISE00";
% Obtain local Earth parameters
GPShLaunchsite = fastinterp2(longitudes, geodeticLatitudes, WGS84ToTerrain, LonLaunch, LatLaunch);
MSLhLaunchsite = fastinterp2(longitudes, geodeticLatitudes, GeoidToTerrain, LonLaunch, LatLaunch);
[ERA0, GMST0, LST0] = getRotAngsfromJDUT1(JDLaunch, deg2rad(LonLaunch));

% Corrections and additional parameters
% 
% SCRIPT (sees entire workspace and adds to/changes it without explicitly
% indicating outputs)
finalizeInputs

% Model the trajectory of a multi-staged rocket launching from the ground
% at some longitude and geodetic latitude as it flies through the
% atmosphere
[t, x, te, xe, ie] = deal(cell(stageNums(end), 1));
for s = stageNums'
    [t{s}, x{s}, te{s}, xe{s}, ie{s}] = ode113(@(tdum, xdum) odeval(tdum, xdum, pars, s), [t0, tf], x0, odeOpts);
    t0 = t{s}(end);
    x0 = x{s}(end, :);
end