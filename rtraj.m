% The (R)ocket (Traj)ectory Program (rtraj)
% Matt Werner (m.werner@vt.edu)

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
Runtime   = true; % true/false - Show runtime after finishing
ODEtime   = true; % true/false - Show solver times at each time step
useCache  = true; % true/false - Utilize a cache to load the environment faster if applicable

%% Inputs
% Obtain actual time right now
tLocalNow = datetime('now', 'timezone', 'America/New_York');
% Obtain actual launch time as Julian date in UT1
[JDLaunch, tLaunchUTC, tLaunch] = defineLaunchTime('15-July-2091 13:53:45.234', 'd-MMMM-yyyy HH:mm:ss.SSS', '-05:00');
% Provide earth model parameters
earthModel = "WGS84"; % Ellipsoid
gravityModel = "tide-free"; % tide-free or zero-tide (EGM08)
gravityDegree = 10; % Degree of gravity model (EGM08)
gravityOrder = 10; % Order of gravity model (EGM08)
magneticDegree = 12; % Degree of magnetic model (WMM20)
magneticOrder = 12; % Order of magnetic model (WMM20)
terrainAngleUnits = "degrees"; % Units for lat/lon on terrain map
atmosModel = "NRLMSISE00"; % Atmosphere model
% Provide launch site location
LonLaunch = convUnits(-119.0560102, "degrees", terrainAngleUnits); % Longitude of launch site [deg or rad]
LatLaunch = convUnits(+ 40.9107330, "degrees", terrainAngleUnits); % Geodetic latitude of launch site [deg or rad]
% Provide launch site temperature and density
localTemp = 70; % Local surface temperature [F]
% Define launch tower length
towerSpan = convUnits(30, "feet", "meters"); % Length of launch tower [m]
% Define rocket parameters (lengths/distances)
numStages = 2; % Number of vehicle stages
LenRocket = [6.35; 3.71]; % Overall length of the rocket during each stage [m]
LenNoseCM = [4.45; 2.71]; % Distance from the nose cone's tip to the center of mass of the fully fueled rocket for each stage [m]stage [m]
diaOuter  = convUnits([8.50; 6.00], "inches", "meters"); % Rocket casing's outermost diameter for each stage [m]
diaThroat = convUnits([2.50; 1.60], "inches", "meters"); % Nozzle throat diameter [m]
diaExit   = convUnits([6.01; 4.47], "inches", "meters"); % Nozzle exit diameter [m]
diaFlatDM = convUnits([36, 108; 36, 108], "inches", "meters"); % Flattened drogue & main parachute diameters (0 indicates no chute) [m]
% Define rocket parameters (moment of inertia)
% The moment of inertia MUST be given relative to the principal axis frame
% centered at the center of mass having its 1-axis point towards the rail
% guides and its 3-axis pointing through the nose
IG(:,:,1) = diag([387.53, 387.54, 1.40]); % Principal moment of inertia
% Define rocket parameters (propulsion)
massInit  = [194.43; 65.33]; % Initial mass of the rocket at the beginning of each burn time [kg]
massMotor = [63.49; 23.58]; % Motor mass of each stage [kg]
burnTimes = [6.95; 5.09]; % Burn time of each stage (has priority over max time in thrust profile) [s]
Yexhaust  = [1.25; 1.25]; % Perfect ratio of specific heats for motor exhaust gases []
FTprofile = ["profiles/prop/Hokie075_Thrust_Stage1.txt", "profiles/prop/S235.csv"];
MFprofile = ["profiles/prop/S1SL.csv", "profiles/prop/S235.csv"];
PCprofile = ["", ""];
% Define rocket parameters (aerodynamics/geometry)
RasAeroCD = {csvread('profiles/drag/H7512.csv', 1, 0); csvread('profiles/drag/H7522.csv', 1, 0)}; % Obtain drag coeff. profile from RasAero csv files
% ... geometry here to define our own CD profile
% Define rocket parameters (staging behavior)
delaySep  = [2; NaN]; % Time it takes for stages to separate after burn-out [s]
delayIgn  = [0; NaN]; % Time it takes for stages to ignite after separation [s]
% Define rocket parameters (controlled descent)
deployAlt = convUnits([2000; 2000], "feet", "meters"); % Altitude of main parachute deployment [m]
% Constants and measured quantities relevant to ODE initial conditions
Rail2Rckt = diaOuter(1)/2; % Perpendicular distance from rail to rocket's centerline [m]
Rail2Vert = 10; % Angle made between railing and vertical [deg]
East2DwnR = 30; % Angle made between launch (downrange) direction and due East [deg]
veps = eps; % Very small initial velocity [m/s]

% SCRIPTS (sees entire workspace and adds to/changes it without explicitly
% indicating outputs)
loadEarthModel % Earth model - gravity & magnetic fields, terrain, and atmosphere
loadPropulsionData
loadCDProfile
finalizeInputs % Corrections and additional parameters

% Model the trajectory of a multi-staged rocket launching from the ground
% at some longitude and geodetic latitude as it flies through the
% atmosphere
[t, x, te, xe, ie] = deal(cell(stageNums(end), 1));
for s = stageNums'
    [t{s}, x{s}, te{s}, xe{s}, ie{s}] = ode113(@(tdum, xdum) odeval(tdum, xdum, pars, s), [t0, tf], x0, odeOpts);
    t0 = t{s}(end);
    x0 = x{s}(end, :);
end