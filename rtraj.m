% OLVT - Trajectory and Analysis
% Simulate ground-to-orbit trajectory of all rocket stages
% Initially created - Oct 16, 2019
% Reboot v7.0.0 - Nov 30, 2020
% Current v7.0.1 - Dec 5, 2020
% Matt Werner (m.werner@vt.edu)

% Notes
% 1. Use a mono font to correctly display formatting of comments/tables
% 2. Disable wrap-lines in Command Window preferences to view tables easily
% 3. No support to integrate automatically on an interval [t0, t1] for t0 ~= 0
% 4. Automatically resets any settings changed by producing a run

% Known problems
%

% To do:
% Everything
% deval() - Can evaluate solution at points (useful for equal time step ?)

% % Suppress warnings induced by switch options
% %#ok<*UNRCH> (unreachable)
% %#ok<*DEFNU> (unused func)
% %#ok<*NASGU> (unused var)
% %#ok<*ALIGN> (misaligned)

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
angleunit = 'rad'; % 'deg'/'rad' - Angle units displayed in plot outputs
USCSIunit = 'SI'; % 'USC'/'SI' - Units in which to display results (US Customary / Metric) (Not yet implemented)
odestats_ = 'on'; % 'on'/'off' - Give ODE solver statistics
oderefine = 3; % Positive integer - Number of points to evaluate between time steps
odereltol = 1e-5; % Positive scalar - Relative error tolerance bound
odeabstol = 1e-7; % Positive scalar - Absolute error tolerance bound
odeinstep = 1e-3; % Positive scalar - Initial step size in time
odemxstep = inf; % Positive scalar - Maximum permissible time step (including inf)
OmitCalcs = false; % true/false - Save time by performing no additional calculations
TinyTable = false; % true/false - Short table showing column names and end of data
FullTable = false; % true/false - Entire table set of all data
DispTable = false; % true/false - Event table
MakePlots = false; % true/false - Show plots
RealTimes = true; % true/false - Show runtime after finishing
ODEsTimes = true; % true/false - Show solver times at each step
HavePrtrb = false; % true/false - Consider perturbations to nominal flight

% Obtain conversions between common units
[ft2in, fps2mph, in2m, lbf2N, lb2kg, psi2Pa] = getUnitConversions;
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
LenRocket = [6.35, 3.71]; % [Overall length during S1, Overall length during S2] [m]
LenNoseCM = 4.15; % Distance from the nose cone's tip to the center of mass [m]
diaOuter_ = [8.5, 6] * in2m; % [S1 outer diameter, S2 outer diameter] [in] --> [m]
diaThroat = [2.5, 1.6] * in2m; % [S1 throat diameter, S2 throat diameter] [in] --> [m]
diaExit__ = [6.01, 4.47] * in2m; % [S1 exit diameter, S2 exit diameter] [in] --> [m]
diaFlatDM = [36, 108] * in2m; % [flat drogue diameter, flat main diameter] [in] --> [m]
% Define rocket parameters (propulsion)
massInit_ = [194.43, 65.33]; % [S1 initial mass, S2 initial mass] [kg]
massMotor = [63.49, 23.58]; % [S1 motor mass, S2 motor mass] [kg]
RSHxhaust = [1.25, 1.25]; % [S1 perfect (R)ratio of (S)pecific (H)eats, S2 perfect ratio of specific heats] []
BurnSimFT = {csvread('csvs/S1SL.csv', 1, 0), csvread('csvs/S235.csv', 1, 0)}; % Obtain data from BurnSim csv files
% Define rocket parameters (aerodynamics)
RasAeroCd = {csvread('csvs/H7512.csv', 1, 0), csvread('csvs/H7522.csv', 1, 0)}; % Obtain data from RasAero csv files
% Define rocket parameters (staging behavior)
delaySep_ = 2; % Time it takes for second stage to separate from first stage [s]
delayIgn_ = 0; % Time it takes for second stage to fire after separation [s]
deployAlt = 609; % Altitude of main parachute deployment (2000 ft) [m]
% Constants and measured quantities relevant to ODE initial conditions
Rail2Rckt = diaOuter_(1)/2; % Perpendicular distance from rail to rocket's centerline [m]
Rail2Vert = 5; % Angle made between railing and vertical [deg]
East2Dwnr = 30; % Angle made between launch (downrange) direction and due East [deg]
veps = 0; % Very small initial velocity [m/s]
% Corrections and additional parameters
% [, , csvs, burnTime, maxburnTime, ...
%     referenceArea, throatArea, exitArea, parachuteArea, hpc] = additional
%  Pick up here

% Obtain parameters that define an Earth model
[GM, Req, Rpo, f, e, w] = defineEllipsoidParameters('WGS84');
Ravg = getEllipsoidAverageRadius(Req, Rpo, '++');
[nG, mG, ~, Cnm, Snm, ~, ~] = loadGravitationalCoefficients(2190, 2159, 'tide-free');
[Cnm, Snm] = updateGravitationalCoefficients(JDLaunch, nG, mG, Cnm, Snm); clear nG mG CnmGo SnmGo
[nM, mM, gnm, hnm, dgnmdt, dhnmdt] = loadMagneticCoefficients(12, 12);
[gnm, hnm] = updateMagneticCoefficients(JDLaunch, gnm, hnm, dgnmdt, dhnmdt); clear dgnmdt dhnmdt
[longitude, geodeticLatitude, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain('radians');

% localElev = fastinterp2(
% ERA0, GMST0, LST0 = getRotAngsfromJDUT1(...)
% stndAtmosTempDiff, air_density0 = ...