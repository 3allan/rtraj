% The (R)ocket (Traj)ectory Program (rtraj)
% Matt Werner (m.werner@vt.edu)

% Notes
% 1. Use a mono font to correctly display formatting of comments/tables
% 2. Disable wrap-lines in Command Window preferences to view tables easily
% 3. Automatically resets any settings changed by producing a run
% 4. Requires Aerospace toolbox

% Stage 1
% 1 - Rail clearance        3 - Separation Delay        5 - Parachute
% 2 - Burnout 1             4 - Apogee                  6 - Touchdown
% Stage 2
% 1 - Rail clearance        3 - Separation Delay        5 - Burnout 2        7 - Parachute
% 2 - Burnout 1             4 - Ignition Delay          6 - Apogee           8 - Touchdown
% ... etc.
% Options
flags.options.display.tablePhases = 1:6; % i:j - Sequences to display in tables (associated with DispTable)
flags.options.display.plotsPhases = 1:6; % i:j - Sequences to display in plots (associated with Plots)
flags.options.ode.stats = 'on'; % 'on'/'off' - Give ODE solver statistics
flags.options.ode.refine = 3; % Positive integer - Number of points to evaluate between time steps
flags.options.ode.relativeTolerance = 1e-5; % Positive scalar - Relative error tolerance bound
flags.options.ode.absoluteTolerance = 1e-7; % Positive scalar - Absolute error tolerance bound
flags.options.ode.initialStep = 1e-3; % Positive scalar - Initial step size in time
flags.options.ode.maxStep = inf; % Positive scalar - Maximum permissible time step (including inf)
% Flags
flags.options.include.perturbations = false; % true/false - Consider perturbations to nominal flight
flags.options.include.cloneDynamics = false; % true/false - Save time by performing no additional calculations
flags.options.show.plots = false; % true/false - Show plots
flags.options.show.runtime = true; % true/false - Show runtime after finishing
flags.options.show.ODEtime = true; % true/false - Show solver times at each time step
flags.options.use.cache = true; % true/false - Utilize a cache to load the environment faster if applicable

%% Inputs
% Obtain actual time right now
time.now = datetime('now', 'timezone', 'America/New_York');
% Obtain actual launch times as Julian date (UT1) and other timezones
% (UTC & local)
[time.launch.JulianDate, time.launch.UTC, time.launch.local] = defineLaunchTime('15-July-2091 13:53:45.234', 'd-MMMM-yyyy HH:mm:ss.SSS', '-05:00');

% Provide earth model parameters
earth.model = "WGS84"; % Ellipsoid
earth.gravity.model = "tide-free"; % tide-free or zero-tide (EGM08)
earth.gravity.degree = 11; % Degree of gravity model (EGM08)
earth.gravity.order = 10; % Order of gravity model (EGM08)
earth.magnetic.degree = 12; % Degree of magnetic model (WMM20)
earth.magnetic.order = 12; % Order of magnetic model (WMM20)
earth.atmos.model = "NRLMSISE00"; % Atmosphere model
earth.terrain.angleUnits = "degrees"; % Units for lat/lon on terrain map
% Provide launch site location
earth.launch.longitude = convUnits(-119.0560102, "degrees", earth.terrain.angleUnits); % Longitude of launch site [deg or rad]
earth.launch.latitude  = convUnits(+ 40.9107330, "degrees", earth.terrain.angleUnits); % Geodetic latitude of launch site [deg or rad]
% Provide launch site temperature and density
launch.temperature = 70; % Local surface temperature [F]
% Define launch tower length
launch.towerSpan = convUnits(30, "feet", "meters"); % Length of launch tower [m]

% Provide a name for the rocket
rocket.name = "Hokie 0.75";
% Define the number of stages in the rocket
rocket.stages = 2; % Number of vehicle stages
% Define rocket parameters (lengths/distances)
rocket.nosecone.series = "Parabolic"; % Series/type of the nosecone (Parabolic/Haack (Haack not yet implemented))
rocket.nosecone.k = 0.9; % Nondimensional parameter defining the nosecone's curvature according to the series
rocket.nosecone.OD = convUnits(6, "inches", "meters"); % Outer diameter of the nosecone (at its base) [m]
rocket.nosecone.length = convUnits(30, "inches", "meters"); % Length of the nosecone from tip to base [m]
rocket.nosecone.thickness = convUnits(0.5, "inches", "meters"); % Distance separating the inside (negative) tip from the outside (positive) tip [m]
rocket.nosecone.density = 2700; % Density of the nosecone material [kg/m3]

rocket.body.OD = convUnits([8.50; 6.00], "inches", "meters"); % Rocket casing's outermost diameter for each stage [m]
rocket.body.ID = convUnits([8.00; 5.50], "inches", "meters"); % Rocket casing's innermost diameter for each stage [m]
rocket.body.length = convUnits([8.2583; 10.37], "feet", "meters"); % Length of each bodytube/casing [m]
rocket.body.density = [2700; 2700]; % Density of the body material [kg/m3]

rocket.frustum.rearOD = convUnits([8.50; 0], "inches", "meters"); % Rocket frustum's outermost diameter for each stage on the face closer to the nozzle [m]
rocket.frustum.frontOD = convUnits([6.00; 0], "inches", "meters"); % Rocket frustum's outermost diameter for each stage on the face closer to the nosecone [m]
rocket.frustum.length = convUnits([6; 0], "inches", "meters"); % Length of each frustum [m]
rocket.frustum.thickness = convUnits([1; 0], "inches", "meters"); % Distance in the axisymmetric direction between the inner and outer faces (NOT the normal distance between the two faces)
rocket.frustum.density = [2700; 0]; % Density of the (conical) frustum material [kg/m3]

rocket.motor.length = [1; 1]; % Length of each motor [m]
rocket.motor.TD = convUnits([2.50; 1.60], "inches", "meters"); % Nozzle throat diameter [m]
rocket.motor.ED = convUnits([6.01; 4.47], "inches", "meters"); % Nozzle exit diameter [m]
rocket.parachute.drogue.diameter = convUnits([36; 36], "inches", "meters"); % Flattened drogue parachute diameter (0 indicates no chute) [m]
rocket.parachute.main.diameter = convUnits([108; 108], "inches", "meters"); % Flattened main parachute diameter (0 indicates no chute) [m]

% Define rocket parameters (propulsion)
profiles.thrust.path = "profiles/prop/" + ["Hokie075_Thrust_Stage1.txt"; "Hokie075_Thrust_Stage2.txt"];
profiles.massFlowRate.path = "profiles/prop/" + ["Hokie075_MassFlowRate_Stage1.txt"; "Hokie075_MassFlowRate_Stage2.txt"];
profiles.burnDepth.path = "profiles/prop/" + ["Hokie075_BurnDepth_Stage1.txt"; "Hokie075_BurnDepth_Stage2.txt"];
profiles.chamberPressure.path = [""; ""];

% Define rocket parameters (aerodynamics/geometry)
% ... geometry here to define our own CD profile

% Define rocket parameters (staging behavior)
rocket.delays.separation = [2; NaN]; % Time it takes for stages to separate after burn-out [s]
rocket.delays.ignition = [0; NaN]; % Time it takes for stages to ignite after separation [s]

% Define rocket parameters (controlled descent)
rocket.parachute.main.deployAtAltitude = convUnits([2000; 2000], "feet", "meters"); % Altitude of main parachute deployment [m]

% Constants and measured quantities relevant to ODE initial conditions
rocket.initial.distance.railToCenterOfMass = rocket.body.OD(1,1)/2; % Perpendicular distance from rail to rocket's centerline [m]
rocket.initial.angle.railToVertical = 10; % Angle made between railing and vertical [deg]
rocket.initial.angle.eastToDownrange = 30; % Angle made between launch (downrange) direction and due East [deg]
rocket.initial.velocity = eps; % Very small initial velocity [m/s]

% SCRIPTS (sees entire workspace and adds to/changes it without explicitly
% indicating outputs)
loadEarthModel % Earth model - gravity & magnetic fields, terrain, and atmosphere
loadPropulsionData % Propulsion model - thrust, mass flow rate & distribution, and (optional) chamber pressure
finalizeInputs % Corrections and additional parameters

%% Flight Model
% Model the trajectory of a multi-staged rocket launching from the ground
% at some longitude and geodetic latitude as it flies through the
% atmosphere
[t, x, te, xe, ie] = deal(cell(stageNums(end), 1));
for s = stageNums'
    [t{s}, x{s}, te{s}, xe{s}, ie{s}] = ode113(@(tdum, xdum) odeval(tdum, xdum, pars, s), [t0, tf], x0, odeOpts);
    t0 = t{s}(end);
    x0 = x{s}(end, :);
end