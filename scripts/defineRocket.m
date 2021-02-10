% 
% Jan. 16, 2021
% 
% Define the flight vehicle as a Matlab structure that holds information
% about its geometry and physical parameters (material/mass density, motor
% characteristics, parachutes, and flight behavior)
% 

% Provide a name for the vehicle
rocket.name = "Hokie 0.75";

% Define the amount of stages in the rocket
rocket.stages = 2; % Number of vehicle stages

%% Staging
% 
% Flight vehicles (rockets) often dispose (stage) casings that hold
% hardware/equipment that isn't needed for the mission any longer to reduce
% the effects of gravity and drag loss on the flight plan. When a vehicle
% is ready to stage, it may choose to wait a period of time before
% permanently disconnecting itself from the empty body components and then
% choose to wait again before firing the next stage. It is possible to
% stage the nosecone from the final sustainer, so there are as many
% options/opportunities to stage as there are "stages" (motors to burn).
% 
% Ex: A two-stage rocket (meaning that it has 1 nosecone initially attached
%     to 2 casings, each holding a motor) that keeps its nosecone attached
%     to the sustainer after the final burnout only stages once. It has the
%     potential to stage once more, but it doesn't by design.
% 
% Ex: A two-stage rocket that stages its nosecone after final burnout
%     stages twice. It does this because the rocket is made up of 3 main
%     components despite being a 2-stage rocket - the two staging bodies
%     and the nosecone.
% 
% A time-to-stage of infinity indicates that the vehicle is finished
% staging; it will never stage past the first instance in which the
% separation delay after burnout is infinity.
% 
% A time-to-ignition of infinity indicates that the vehicle is finished
% firing; it will never fire another (main) motor past the first instance
% in which the ignition delay (after performing a separation) is infinity.
% 

% Time for stages to separate after burnout [s]
rocket.delays.separation(2,1) = inf;
rocket.delays.separation(1,1) = 2;

% Time for stages to ignite after separation [s]
rocket.delays.ignition(2,1) = inf;
rocket.delays.ignition(1,1) = 2;

%% Parachutes
% 
% Note: Each stage usually contains at least 2 parachutes (not counting
% backups) - the drogue parachute and the main parachute. The drogue
% parachute usually deploys at apogee, though may also be deployed at some
% target altitude (always while falling back down to Earth's surface), to
% have the effect of keeping the vehicle upright during decent. Such a
% descent is confirmed by checking that  the vehicle's angle of attack is
% roughly 180 degrees during descent.
% 
% Each parachute has a maximum diameter when laid and spread out on a flat
% surface. This diameter is related to the diameter of a "fully-inflated"
% hemisphere by equating surface areas of a flat disk and hemisphere, which
% is how the parachute is modelled to be when fully deployed. The projected
% area onto the plane that's normal to the velocity direction (i.e. how the
% parachute opens in time) is specified as an explicit function of time and
% NOT dynamically solved.
% 

% Flattened drogue parachute diameter (0 indicates no chute) [m]
rocket.parachute.drogue.flatDiameter(2,1) = convUnits(36, "inches", "meters");
rocket.parachute.drogue.flatDiameter(1,1) = convUnits(36, "inches", "meters");

% Flattened main parachute diameter (0 indicates no chute) [m]
rocket.parachute.main.flatDiameter(2,1) = convUnits(108, "inches", "meters");
rocket.parachute.main.flatDiameter(1,1) = convUnits(108, "inches", "meters");

% Parachute-opening profile. This profile determines the percentage of how
% much of the parachute has inflated since deployment (i.e. how much area
% it is projecting onto the plane normal to the velocity vector since it
% was ejected and is opening in time). As such, the profile should
% MONOTONICALLY vary from 0 to 1 for all time t > 0. Real parachutes take
% some time trailing behind before beginning to open, and the actual
% opening motion is fast, so a realistic transition from closed (0) to
% fully open (1) in terms of a mathematical function means that the
% functions perhaps stays near 0 for some time followed by a sharp gradient
% to 1 where it permanently stays.
% 
% The nosecone gets its own option for a parachute in case the nosecone
% detaches from the final sustainer as part of the flight plan. Therefore,
% there is 1 more entry than the number of stages for the parachute
% profiles.
% 
% Ex: A 2-stage rocket that doesn't stage its nosecone has a profile for
% its booster, sustainer, and its nosecone. The nosecone 
% 
% No checks are performed to determine if the function varies monotonically
% from 0 to 1. (Therefore, it's possible to define a parachute that
% repeatedly opens and closes in time during its descent.)
% 
% Drogue parachute profiles
rocket.parachute.drogue.openingFunc{3,1} = @(t) 0;
rocket.parachute.drogue.openingFunc{2,1} = @(t) erf(5*t);
rocket.parachute.drogue.openingFunc{1,1} = @(t) erf(5*t);
% Main parachute profiles
rocket.parachute.main.openingFunc{3,1} = @(t) 0;
rocket.parachute.main.openingFunc{2,1} = @(t) heaviside(t-2).*erf(2*t-4);
rocket.parachute.main.openingFunc{1,1} = @(t) heaviside(t-2).*erf(2*t-4);

% Indication as to how to interpret the distance/input at which the drogue
% and main parachutes will be deploying
% Available options are:
%  - "ASL" or "Above Sea Level"
%  - "AGL" or "Above Ground Level"
%  - "APO" or "Apogee"
%  - "BAP" or "Below Apogee"
%  - "ATP" or "Atmospheric Pause"
% 
% Drogue parachute
rocket.parachute.drogue.deploymentType(3, 1) = "Apogee";
rocket.parachute.drogue.deploymentType(2, 1) = "Apogee";
rocket.parachute.drogue.deploymentType(1, 1) = "Apogee";
% Main parachute
rocket.parachute.main.deploymentType(3, 1) = "AGL";
rocket.parachute.main.deploymentType(2, 1) = "AGL";
rocket.parachute.main.deploymentType(1, 1) = "AGL";

% Altitude of drogue parachute deployment [m]
% Permissible options are:
%  -    "Above Sea Level": z > 0
%  - "Above Ground Level": z > local elevation of area
%  -             "Apogee": N/A
%  -       "Below Apogee": tallest mountain < z < apogee
%  -  "Atmospheric Pause": N/A
% 
% Drogue parachute
rocket.parachute.drogue.deploymentAlt(3, 1) = NaN;
rocket.parachute.drogue.deploymentAlt(2, 1) = NaN;
rocket.parachute.drogue.deploymentAlt(1, 1) = NaN;
% Main parachute
rocket.parachute.main.deploymentAlt(3, 1) = 610;
rocket.parachute.main.deploymentAlt(2, 1) = 610;
rocket.parachute.main.deploymentAlt(1, 1) = 610;

%% Nosecone
% 
% Specify the precise shape of the rocket's nosecone using some predefined
% shapes and models. 
% 
% Shaping the nosecone is important for structural AND aerodynamic
% propertiies as a nosecone's profile has a large influence on the rocket's
% drag profile (especially at supersonic speeds due to how shockwaves
% interact with the atmospheric conditions and nosecone geometry).
% 

% Series/type of the nosecone 
% Available options are:
%  - "Linear"
%  - "Power"
%  - "Parabolic"
%  - "Tangent-Ogive"
%  - "Secant-Ogive"
%  - "Haack"
rocket.nosecone.series = "Haack"; 

% Interpretation of thickness
% Available options are:
%  - "Frontal"
%  - "Vertical"
%  - "Vetical-Pointwise"
%  - "Normal" (Not yet implemented)
rocket.nosecone.thicknessType = "Vertical";

% Cone constant (only takes effect for parameterizable series)
% Permissible options are:
%  -        "Linear": N/A
%  -         "Power": 0 <  k <= 1
%  -     "Parabolic": 0 <= k <= 1
%  - "Tangent-Ogive": N/A
%  -  "Secant-Ogive": k > 0.5*(R^2 + L^2)/R
%  -         "Haack": 0 <= k <~ 2/3 (no upper bound but 2/3 is tangent)
rocket.nosecone.k = 0;

% Outer diameter of the nosecone at its base [m]
rocket.nosecone.OD = convUnits(6, "inches", "meters"); 

% Length of the nosecone from tip to base [m]
rocket.nosecone.length = convUnits(30, "inches", "meters");

% Distance separating the inner/outer surfaces according to the specified
% thickness type [m]
% Note: - "Vertical-Pointwise" thickness expects an n-by-2 array specifying
%         the thickness at various points along the nosecone (from tip to
%         base). The first column must be positions along the nosecone in
%         ascending order and the second column must be the associated
%         thicknesses at those respective positions.
%       - Otherwise, the thickness is a single number.
rocket.nosecone.thickness = convUnits(0.25, "inches", "meters");

% Material density of the nosecone [kg/m3]
% Note: The density may be given either as a scalar or as an n-by-2 array.
% ---   Scalar: The density is assumed constant over the entire nosecone.
%   ###               Ex: MDF = 2700;
%   ###               Ex: MDF = defineMassDensityFunction(2700);
%   ###               Ex:   L = rocket.nosecone.length;
%                         MDF = wrapMDF(@(x) 2700, L, 1, "Linear");
% 
%               Note: All 3 of these examples return the exact same
%               quantity - they all return 'MDF = 2700'. The functions in
%               the latter 2 examples are designed to handle varying values
%               of density, so using them to define a constant scalar
%               density is cumbersome.
% 
% ---    Array: The density is provided at various points along the
%               nosecone (from tip to base). The first column must be the
%               positions along the nosecone in ascending order and the
%               second column must be the associated material density at
%               those respective positions. This specification of the
%               density is the "mass density function" (MDF).
%   ###               Ex: L = rocket.nosecone.length;
%                         f = @(x) 2700*cos(x/L/2);
%                         x = linspace(0, L, 200)';
%                       MDF = [x, f(x)];
%               Note: The mass density function may be specified as a
%                     constant by simply providing 2 points with the
%                     positions equivalent to the tip (0) and the base (L)
%                     while the density remains a constant value.
%   ###               Ex:   L = rocket.nosecone.length;
%                           x = [0; L];
%                           d = 2700*[1; 1];
%                         MDF = [x, d];
%               Note: The density may be provided with uncertainty directly
%                     added into it (a generalization of perfect values
%                     with uncertainty equivalent to zero). This
%                     uncertainty can represent a change in material,
%                     cargo/payloads, etc.
%   ###               Ex: n = 200; % Resolution
%                         L = rocket.nosecone.length;
%                         x = linspace(0, L, n)';
%                         u = 0.10; % Relative uncertainty (10%)
%                         U = u*(2*rand(n,1) - 1)
%                         d = 2700*(1 + U);
%                       MDF = [x, d];
%               Note: More customization may be added (including 'perfect'
%                     (no uncertainty) discontinuities) by running the mass
%                     density profile through the following function.
%   ###               Ex: L = rocket.nosecone.length;
%                         x = linspace(0, L, 500)';
%                         f = 2700*cos(x/L/2);
%                         u = [0.01; 0.075; 0.1; 0.1];
%                         I = [0, 1; 0.4, 0.5; 0.7, 0.8; -0.01, 0.01];
%                         X = [0, L; 0, L/30; L/15, L/3.75; L/1.5, L];
%                       MDF = defineMassDensityFunction(f, x, u, I, X);
%               Note: The workspace can remain clean of temporary variables
%                     by calling this previous function using the wrapper.
%   ###               Ex: MDF = wrapMDF(@(x) 2700 + x, ...
%                                 rocket.nosecone.length, ...
%                                 250, ...
%                                 "Linear", ...
%                                 [0.3; 0.95], ...
%                                 [-1, 1; 0.8, 0.8], ...
%                                 [0, 0.2; 0.3,0.6]*rocket.nosecone.length)
% 
% (rocket.nosecone.density = MDF)
% 
% Density function of the nosecone [m | kg/m3]
rocket.nosecone.density = ...
    wrapMDF(@(x) 2700*cos(x/rocket.nosecone.length/2), ...
            rocket.nosecone.length, ...
            500, ...
            "Cosine", ...
            [0.03; 1.075; 0.1275; 0.1], ...
            [0, 1; 0.5, 0.5; 1, 1; -0.01, 0.01], ...
            [0, 1; 0, 1/30; 1/10, 1/1.2; 1/1.5, 1]);

%% Shoulders (also called connectors, frustums, etc.)
% 
% Note: The nosecone is permitted to stage away from the final sustainer if
% the nosecone is connected to the final sustainer by a shoulder (no matter
% how small). 
% 
% To keep the nosecone and final sustainer permanently attached after the
% final burnout, set all of the parameters for the shoulder connecting the
% nosecone to the final sustainer to 0. This shoulder exists in the
% position equal to the number of stages for the vehicle.
% 
% Ex: A 2-stage rocket that keeps its nosecone attached to the body after
%     the final burnout only has 1 shoulder to stage the booster.
%     Therefore, such a rocket doesn't have an interface for separating the
%     nosecone from the sustainer, so the shoulder parameters (in position
%     2) are all 0.
% 
% Ex: A 2-stage rocket that separates its nosecone from the sustainer after
%     the final burnout has 2 shoulders - 1 to stage the booster and 1 to
%     stage the sustainer. In this case, the shoulder parameters (in
%     position 2) should reflect a physical shoulder piece (i.e.,
%     everything is not 0). The shoulder IS permitted to be a straight
%     cylinder, so a nosecone fitted into the final sustainer could be
%     modelled as a slightly shortened casing following by a cylindrical
%     shoulder.
% 
% Note: These pieces are meant only to connect body components in an
% increasing-diameter fashion from the nosecone downwards. In other words,
% these components cannot currently be used as boattails.
% 

% Thickness type specification (Vertical/Normal)
rocket.shoulder.thicknessType(2,1) = "Vertical";
rocket.shoulder.thicknessType(1,1) = "Normal";

% Outer diameter of the shoulder at its smaller end (facing towards the nosecone) [m]
rocket.shoulder.frontOD(2,1) = 0;
rocket.shoulder.frontOD(1,1) = convUnits(6.00, "inches", "meters");

% Outer diameter of the shoulder at its larger end (facing away from the nosecone) [m]
rocket.shoulder.rearOD(2,1) = 0;
rocket.shoulder.rearOD(1,1) = convUnits(8.50, "inches", "meters");

% Distance between the shoulder's two faces [m]
rocket.shoulder.length(2,1) = 0;
rocket.shoulder.length(1,1) = convUnits(6.00, "inches", "meters");

% Distance separating the inner/outer surfaces according to the specified
% thickness type [m]
rocket.shoulder.thickness(2,1) = 0;
rocket.shoulder.thickness(1,1) = convUnits(0.50, "inches", "meters");

% Material density of the shoulder(s) [kg/m3]
% Note: The density may be given either as a scalar or as an n-by-2 array
%       for each shoulder piece. The information, regardless of being a
%       scalar or an array, needs to be stored in a cell array per stage.
% 
% ---   Scalar: The density is assumed constant over this stage's frustum.
%   ###             Ex: MDF = 2700;
%   ###             Ex: MDF = defineMassDensityFunction(2700);
% 
% ---    Array: The density is provided at various points along the frustum
%               (from the smaller OD face to the larger OD face). The first
%               column must be the positions along the nosecone in
%               ascending order and the second column must be the
%               associated material density at those respective positions.
%               This specification of the density is the "mass density
%               function" (MDF).
%   ###               Ex: L = rocket.shoulder.length(1,1);
%                         f = @(x) 2700*(1 + exp(-(x - 0.3*L)^2));
%                         x = linspace(0, L, 200)';
%                       MDF = [x, f(x)];
%               Note: The density may be provided with uncertainty directly
%                     added into it (a generalization of perfect values
%                     with uncertainty equivalent to zero). This
%                     uncertainty can represent a change in material,
%                     cargo/payloads, etc.
%   ###               Ex: L = rocket.shoulder.length(1,1);
%                         x = linspace(0, L, 500)';
%                         f = 2700*(1 + exp(-(x - 0.3*L).^2));
%                         u = [0.05; 0.075];
%                         I = [-0.3, 0.1; 0.7, 0.8];
%                         X = [L/30; L/25, L/3.75; L/1.5];
%                       MDF = defineMassDensityFunction(f, x, u, I, X);
% 
% See the nosecone mass density function description for more examples. The
% only difference between the implementation of this MDF and the nosecone's
% MDF is that the results must go into a cell array corresponding to the
% frustum's stage.
% 
% (rocket.shoulder.density{stage, 1} = MDF;)
% 
% Density function of the nosecone material [m | kg/m3]
rocket.shoulder.density{2,1} = 0;
rocket.shoulder.density{1,1} = wrapMDF(@(x) 2700*(1 + exp(-x.^2)), ...
                                       rocket.shoulder.length(1,1), ...
                                       50, ...
                                       "Linear", ...
                                       [0.05; 0.075], ...
                                       [-0.3, 0.1; 0.7, 0.8], ...
                                       [1/30, 1/25; 1/3.75, 1/1.5]);

%% Cylinders (also called casings, body/tubes, etc.)
% 
% Note: These components are simply hollow cylinders; curvature in the
% outermost body is not permitted. The cylinder must be able to fit the
% solid motor inside.
% 

% Rocket casing's outermost diameter for each stage [m]
rocket.cylinder.OD(2,1) = convUnits(6.00, "inches", "meters"); 
rocket.cylinder.OD(1,1) = convUnits(8.50, "inches", "meters"); 

% Rocket casing's innermost diameter for each stage [m]
rocket.cylinder.ID(2,1) = convUnits(5.50, "inches", "meters"); 
rocket.cylinder.ID(1,1) = convUnits(8.00, "inches", "meters"); 

% Length of each bodytube/casing [m]
rocket.cylinder.length(2,1) = convUnits(10.37, "feet", "meters"); 
rocket.cylinder.length(1,1) = convUnits(8.2583, "feet", "meters"); 

% Material density of the casing(s) [kg/m3]
% Note: The density may be given either as a scalar or as an n-by-2 array
%       for each shoulder piece. The information, regardless of being a
%       scalar or an array, needs to be stored in a cell array per stage.
% 
% See the frustum MDF description for help.
rocket.cylinder.density{2,1} = wrapMDF(@(x) 2700, ...
                                       rocket.cylinder.length(2,1), ...
                                       50, ...
                                       "Linear", ...
                                       0.2, ...
                                       [1, 1], ...
                                       [0.1, 0.4]);
rocket.cylinder.density{1,1} = wrapMDF(@(x) 2700, ...
                                       rocket.cylinder.length(1,1), ...
                                       50, ...
                                       "Linear", ...
                                       0.2, ...
                                       [1, 1], ...
                                       [0.1, 0.4]);
                                   
%% Fins
% Fin geometry here (to define CD profile, provide spinning moments, etc.)

%% Motors
% 
% Motor characteristics as functions of burn time are provided by custom
% csv-style .txt files. Additional information regarding the profiles,
% particularly thrust, are also found in these text files along with the
% profiles. Such additional information within the thrust profile includes
% the ambient pressure at which the thrust curve was generated.
% 

% Path(s) to the motor thrust profile(s)
rocket.motor.files.thrust(2,1) = "profiles/prop/Hokie075_Thrust_Stage2.txt";
rocket.motor.files.thrust(1,1) = "profiles/prop/Hokie075_Thrust_Stage1.txt";

% Path(s) to the motor mass flow rate profile(s)
rocket.motor.files.massFlowRate(2,1) = "profiles/prop/Hokie075_MassFlowRate_Stage2.txt";
rocket.motor.files.massFlowRate(1,1) = "profiles/prop/Hokie075_MassFlowRate_Stage1.txt";

% Path(s) to the motor burn depth profile(s)
rocket.motor.files.burnDepth(2,1) = "profiles/prop/Hokie075_BurnDepth_Stage2.txt";
rocket.motor.files.burnDepth(1,1) = "profiles/prop/Hokie075_BurnDepth_Stage1.txt";

% Path(s) to the motor chamber pressure profile(s)
rocket.motor.files.chamberPressure(2,1) = "";
rocket.motor.files.chamberPressure(1,1) = "";

% Mass [kg]
rocket.motor.mass(2,1) = 23.58;
rocket.motor.mass(1,1) = 63.49;

% Length [m]
rocket.motor.length(2,1) = convUnits(49.25, "inches", "meters");
rocket.motor.length(1,1) = convUnits(60.00, "inches", "meters");

% Outer-diameter [m]
rocket.motor.OD = rocket.cylinder.ID;

% Inner-diameter [m]
rocket.motor.ID(2,1) = convUnits(1.60, "inches", "meters");
rocket.motor.ID(1,1) = convUnits(2.50, "inches", "meters");

%% Nozzles
% Nozzle geometry here (adds more weight, specifies nozzle geometry
% (throat, exit), direction of thrust (perturbations), etc.)

% Nozzle mass [kg]
rocket.nozzle.mass(2,1) = convUnits( 8.09, "lb", "kg");
rocket.nozzle.mass(1,1) = convUnits(21.30, "lb", "kg");

% Nozzle length [m]
rocket.nozzle.length(2,1) = convUnits(7.83, "inches", "meters");
rocket.nozzle.length(1,1) = convUnits(9.65, "inches", "meters");

% Nozzle throat diameter [m]
rocket.nozzle.TD(2,1) = convUnits(1.60, "inches", "meters"); 
rocket.nozzle.TD(1,1) = convUnits(2.50, "inches", "meters"); 

% Nozzle exit diameter [m]
rocket.nozzle.ED(2,1) = convUnits(4.47, "inches", "meters");
rocket.nozzle.ED(1,1) = convUnits(6.01, "inches", "meters");

%% Initial configuration
% 
% Constants and measurable quantities of the vehicle while set up on the
% launch rail. These quantities are relevant to ODE's initial conditions.
% 

% Perpendicular distance from rail to rocket's centerline [m]
rocket.initial.railOffset = 0.05 + max(rocket.cylinder.OD)/2;

% Angle made between railing and vertical [deg]
rocket.initial.railToVertical = 10;

% Angle made between launch (downrange) direction and due East [deg]
rocket.initial.eastToDownrange = 30;

% Very small initial velocity [m/s]
rocket.initial.speed = 0; 