% OLVT - Trajectory and Analysis
% Simulate ground-to-orbit trajectory of all rocket stages
% Initially created - Oct 16, 2019
% Reboot v7.0.0 - Nov 30, 2020
% Current v7.0.7 - Dec 11, 2020
% Matt Werner (m.werner@vt.edu)
clearvars

% Notes
% 1. Use a mono font to correctly display formatting of comments/tables
% 2. Disable wrap-lines in Command Window preferences to view tables easily
% 3. Automatically resets any settings changed by producing a run
% 4. Requires Aerospace toolbox

% Known problems
% - Need to finish first before evaluating problems

% To do:
% 1. Consider fully replacing fastinterp2 and the terrain grids with simply
%    the harmonic coefficients
% 
% 2. Evaluate the dynamics by finding the various quantities required to
%    evaluate f(t, x) (like aoa, CD, L, torques, etc.) and precompute their
%    expressions (on paper) to increase efficiency
% 
% 3. deval() - Can evaluate solution at points (useful for equal time step ?)

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
tLocalNow = datetime(datetime, 'timezone', 'America/New_York');
% Obtain actual launch time as Julian date in UT1
[JDLaunch, tLaunchUTC, tLaunch] = defineLaunchTime('15-July-2091 13:53:45.234', 'd-MMMM-yyyy HH:mm:ss.SSS', '-05:00');
% Provide launch site location
LonLaunch = -119.0560102; % Longitude of launch site [deg]
LatLaunch = + 40.910733; % Geodetic latitude of launch site [deg]
% Provide launch site temperature
localTemp = 70; % Local surface temperature [F]
% Define launch tower length
towerSpan = 30 * ft2in * in2m; % Length of launch tower [ft] --> [in] --> [m]
% Define rocket parameters (lengths/distances)
stageNums = (1:2)'; % Number of vehicle stages
LenRocket = [6.35; 3.71]; % [Overall length during S1, Overall length during S2] [m]
LenNoseCM = [4.15; 1.5]; % Distance from the nose cone's tip to the center of mass [m]
diaOuter_ = [8.5; 6] * in2m; % [S1 outer diameter, S2 outer diameter] [in] --> [m]
diaThroat = [2.5; 1.6] * in2m; % [S1 throat diameter, S2 throat diameter] [in] --> [m]
diaExit__ = [6.01; 4.47] * in2m; % [S1 exit diameter, S2 exit diameter] [in] --> [m]
diaFlatDM = [36, 108; 36, 108] * in2m; % [flat drogue diameter, flat main diameter] [in] --> [m]
% Define rocket parameters (propulsion)
massInit_ = [194.43; 65.33]; % [S1 initial mass, S2 initial mass] [kg]
massMotor = [63.49; 23.58]; % [S1 motor mass, S2 motor mass] [kg]
burnTimes = [6.95; 5.09]; % [S1 burn time, S2 burn time] (takes precedence over max time in FTProfile) [s]
FTProfile = {stripcsv('csvs/S1SL.csv', [1, 4], 1); stripcsv('csvs/S235.csv', [1, 4], 1)}; % Obtain data from BurnSim csv files
FTUnits__ = {["s", "lbf"]; ["s", "lbf"]}; % Units
MFProfile = {stripcsv('csvs/S1SL.csv', [1, 6], 1); stripcsv('csvs/S235.csv', [1, 6], 1)}; % Mass flow rate
MFUnits__ = {["s", "lb/s"]; ["s", "lb/s"]}; % Units
PCProfile = {NaN; NaN}; % Chamber pressure 
PCUnits__ = {["NaN", "NaN"]; ["NaN", "NaN"]}; % Units
Yexhaust_ = [1.25; 1.25]; % [S1 perfect ratio of specific heats, S2 perfect ratio of specific heats] []
% Define rocket parameters (aerodynamics/geometry)
RasAeroCd = {csvread('csvs/H7512.csv', 1, 0); csvread('csvs/H7522.csv', 1, 0)}; % Obtain data from RasAero csv files
% ... geometry here
% Define rocket parameters (staging behavior)
delaySep_ = [2; NaN]; % Time it takes for second stage to separate from first stage [s]
delayIgn_ = [0; NaN]; % Time it takes for second stage to fire after separation [s]
% Define rocket parameters (controlled descent)
deployAlt = [609; 609]; % Altitude of main parachute deployment (2000 ft) [m]
% Constants and measured quantities relevant to ODE initial conditions
Rail2Rckt = diaOuter_(1)/2; % Perpendicular distance from rail to rocket's centerline [m]
Rail2Vert = 5; % Angle made between railing and vertical [deg]
East2DwnR = 30; % Angle made between launch (downrange) direction and due East [deg]
veps = 0; % Very small initial velocity [m/s]

% Obtain parameters that define an Earth model
[GM, Req, Rpo, f, e, w] = defineEllipsoidParameters('WGS84');
[nG, mG, ~, Cnm, Snm, ~, ~] = loadGravitationalCoefficients(10, 5, 'tide-free');
[Cnm, Snm] = updateGravitationalCoefficients(JDLaunch, nG, mG, Cnm, Snm); clear CnmGo SnmGo
[nM, mM, gnm, hnm, dgnmdt, dhnmdt] = loadMagneticCoefficients(12, 12);
[gnm, hnm] = updateMagneticCoefficients(JDLaunch, gnm, hnm, dgnmdt, dhnmdt); clear dgnmdt dhnmdt
[longitudes, geodeticLatitudes, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain('degrees');
atmosModel = "NRLMSISE-00";
% Obtain local Earth parameters
ElevationGPS = fastinterp2(longitudes, geodeticLatitudes, WGS84ToTerrain, LonLaunch, LatLaunch);
ElevationMSL = fastinterp2(longitudes, geodeticLatitudes, GeoidToTerrain, LonLaunch, LatLaunch);
[ERA0, GMST0, LST0] = getRotAngsfromJDUT1(JDLaunch, deg2rad(LonLaunch));

% Corrections and additional parameters
% 
% SCRIPT (sees entire workspace and adds to/changes it without explicitly
% indicating outputs)
correctionsToInputs

% Define ODE initial conditions while rocket is resting on the launch rail
% before motor ignition
[x0, v0, q0, w0] = deal(zeros(3, 1));
q0(4, 1) = 1;
xx0 = [x0; v0; q0; w0];
% Define the initial time at which integration begins
t0 = 0;
% Define the final time at which integration ends
tf = 2000;

% Place necessary variables into pars categorized by application and amount
% of rows (tables require equal amounts of rows per structure) (structures
% are accessed by . notation (Ex: pars.time.tLocalNow.TimeZone indicates
% that 'pars' has a structure called 'time' which holds a value called
% 'tLocalNow' which has a property called 'TimeZone')).
% --- ODE ---
pars.options = table(HavePrtrb, CloneData, ShowPlots, Runtime__, ODEtime__);
% 
% --- TIME ---
pars.time = table(tLocalNow, JDLaunch, tLaunchUTC, tLaunch);
% 
% --- LAUNCH ---
% Pass quantities to do with the launch time and initial orientation of the
% rail with respect to the ENV frame.
pars.launchsite = table(LonLaunch, LatLaunch, LonLaunchrad, LatLaunchrad, ... 
                        ERA0, GMST0, LST0, ElevationGPS, ElevationMSL, ...
                        localTemp, towerSpan, w_ecef_rowvec, w_env_rowvec, ...
                        rLaunch_ecef_rowvec);
% 
pars.launchrail = table(Rail2Rckt, Rail2Vert, East2DwnR, veps);
%
% --- FLIGHT VEHICLE ---
% (geometric/phase-defining) together
pars.rocket = table(stageNums, LenRocket, LenNoseCM, ...
                    diaOuter_, diaThroat, diaExit__, diaFlatDM, ...
                    massInit_, massMotor, Yexhaust_, burnTimes, ...
                    delaySep_, delayIgn_, deployAlt);
% 
% --- PROPULSION ---
% thrust
pars.thrustProfile = table(FTProfile, FThpc);
% 
% mass flow rate
pars.massFlowRate = table(MFProfile, MFhpc);
% 
% chamber pressure
pars.chamberPressure = table(PCProfile, PChpc);
% 
% --- AERODYNAMICS ---
% drag coefficient
pars.dragCoefficient = table(RasAeroCd);
% 
% --- EARTH MODEL ---
% ellipsoid parameters
pars.ellipsoid = table(GM, Req, Rpo, f, e, w);
% 
% longitudes of height maps
pars.terrain.longitudes = table(longitudes);
% 
% geodetic latitudes of height maps
pars.terrain.geodeticLatitudes = table(geodeticLatitudes);
% 
% height maps
pars.terrain.height = table(WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain);
% 
% gravitational field potential harmonic coefficients
pars.coefficients.gravity = table(nG, mG, Cnm, Snm);
% 
% magnetic field potential harmonic coefficients
pars.coefficients.magnetf = table(nM, mM, gnm, hnm);
% 
% atmosphere model
pars.atmosphericModel = table(atmosModel);
%
% rotations
pars.rotations = table(Tenv_ecf, Tecf_env);
% -----------------------------------------------------------------------

% Set options for ODE solver
options = odeset('reltol', odereltol, 'abstol', odeabstol, ...
                 'initialstep', odeinstep, 'maxstep', odemxstep, ...
                 'stats', odestats_, 'refine', oderefine, ...
                 'events', @odevents);

% Model the trajectory of a multi-staged rocket launching from the ground
% at some longitude and geodetic latitude as it flies through the
% atmosphere
[t, x, te, xe, ie] = deal(cell(stageNums(end), 1));
for s = stageNums
    [t{s}, x{s}, te{s}, xe{s}, ie{s}] = ode113(@(tdum, xdum) odeval(tdum, xdum, pars), [t0, tf], xx0);
    t0 = t{s}(end);
    xx0 = x{s}(end, :);
end