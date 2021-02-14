% One-time run script that calculates some additional parameters from the
% inputs of rtraj

%% Burn time
% Explicitly place the burn times of the motors in with the thrust
% profile, where the first column indicates the nominal burn time and
% the second column indicates the (absolute) uncertainty
% (This fits better in with 'Propulsion', but it's needed to evaluate a
% term used in computing the wet properties.)
for stage = 1:rocket.stages
    rocket.motor.burnTime(stage, 1:2) = ...
        [rocket.motor.profiles.thrust.curve{stage, 1}(end, 1), 0];
end

%% Structures

% Start the burn depth profile at t = 0
rocket.motor.IDatt = rocket.motor.ID;

rocket = computeDryProperties(rocket);
rocket = computeWetProperties(rocket);


%% Propulsion

% Allocate memory to determining the retarded (as in delayed, but computed
% rather than chosen) times to ignition relative to launch (t = 0)
rocket.motor.retardedTimes = NaN(rocket.stages, 2);

% Objectives:
%  1. Determine delays to firing from launch
%  2. Provide a (linear) mass flow rate if a mass flow rate profile wasn't
%     explicitly given
%  2a. Ensure that the mass flow rate is negative
%  3. Calculate the mass profile
%  4. Provide a (linear variation of) mass distribution if a burn depth
%     profile wasn't explicitly given
%  5. Indicate whether a chamber pressure profile was provided or not
%  6. Determine the 'pchip' coefficients for consistently fast
%     interpolation
% 
% Note: Units for propulsive characteristics are already in SI, so no unit
%       conversions need to be done, and burn times are assumed to be when
%       the thrust reaches 0.
for stage = 1:rocket.stages
    
    % ======================= Retarded Time (1) ===========================
    if (stage == 1)
        % No delays during ignition since ignition occurs at t = 0
        rocket.motor.retardedTimes(1, :) = zeros(1, 2);
    else
        % Provide a nonzero offset due to burntimes, separation delays, and
        % ignition delays of previous stages
        rocket.motor.retardedTimes(stage, 1) = ...
            rocket.motor.retardedTimes(stage-1, 1) + ...
            rocket.motor.burnTime(stage-1, 1) + ...
            rocket.delays.separation(stage-1, 1) + ...
            rocket.delays.ignition(stage-1, 1);
    end
    
    
    % =================== Mass Flow Rate Profile (2) ======================
    % Ensure that the mass flow rate profile is negative (3a)
    if (flags.exists.profile.massFlowRate(stage, 1))
        % Negate the provided mass flow rate profile 
        rocket.motor.profiles.massFlowRate.curve{stage, 1}(:, 2) = ...
            -abs(rocket.motor.profiles.massFlowRate.curve{stage, 1}(:, 2));
    else
        % Define temporary variables for brevity and clarity
        tmp_burnTime = rocket.motor.burnTime(stage, 1);
        tmp_motorMass = rocket.motor.mass(stage, 1);
        tmp_linearMassFlowRate = -abs(tmp_motorMass/tmp_burnTime);
        tmp_ltb = [0; 1]*tmp_burnTime;
        tmp_lmd = [1; 1]*tmp_linearMassFlowRate;
        % Define the linear mass flow rate profile
        rocket.motor.profiles.massFlowRate.curve{stage, :} = [tmp_ltb, tmp_lmd];
        clear tmp_*
    end
    
    
    % ======================== Mass Profile (3) ===========================
    % Calculate the motor's mass as a function of burntime explicitly from
    % the mass flow rate via numerical integration
    
    % Define temporary variables for brevity and clarity
    tmp_t = rocket.motor.profiles.massFlowRate.curve{stage, 1}(:, 1);
    tmp_f = rocket.motor.profiles.massFlowRate.curve{stage, 1}(:, 2);
    
    % Create mass profile from mass flow rate profile and prepare it to be
    % pchipped
    tmp_expelledMass = -cumtrapz(tmp_t, tmp_f);
    
    % Check that the final amount doesn't exceed the specified motor mass
    tmp_excessMass = tmp_expelledMass(end) - rocket.motor.mass(stage,1);
    if (tmp_excessMass > 1.000)
        error("Integrated motor mass exceeds indicated motor mass by %1.2f kg (stage %1.0f)", tmp_excessMass, stage)
    elseif (tmp_excessMass > 0.050)
        warning("Integrated motor mass exceeds indicated motor mass by %1.2f kg (stage %1.0f)", tmp_excessMass, stage)
    elseif (tmp_excessMass < -0.050)
        warning("Indicated motor mass exceeds integrated motor mass by %1.2f kg (stage %1.0f)", -tmp_excessMass, stage)
    elseif (tmp_excessMass < -1.000)
        error("Indicated motor mass exceeds integrated motor mass by %1.2f kg (stage %1.0f)", -tmp_excessMass, stage)
    end
    
    % Define the mass profile
    rocket.mass{stage, 1} = [tmp_t, rocket.wet.mass(stage) - tmp_expelledMass];
    
    
    % ===================== Burn Depth Profile (4) ========================
    % Nothing to do if the profile is provided, but if it's not provided,
    % then make a linear variation varying from the rocket's centerline to
    % its outer diameter0
    if (~flags.exists.profile.burnDepth(stage, 1))
        % Temporary variable for brevity and clarity
        tmp_ID = rocket.cylinder.ID(stage, 1);
        tmp_tb = [0; 1]*rocket.motor.burnTime(stage, 1);
        tmp_bd = [0; 1]*tmp_ID;
        profiles.burnDepth.func{stage, 1} = [tmp_tb, tmp_bd];
    end
    
    % ================== Chamber Pressure Profile (5) =====================
    % Nothing to do since the flags regarding existence of a chamber
    % pressure profile have already been determined; this objective serves
    % only as a reminder that it has already been done.
    
    
    % Clean up temporary variables
    clear tmp_*
    
    
    % =====================================================================
    % = Piecewise Cubic Hermite Interpolating Polynomial Coefficients (6) =
    % =====================================================================
    % Obtain (H)ermite (P)olynomial (C)oefficients (HPC) for fast
    % interpolation ('pchip') of the rocket's mass, thrust, burn depth, and
    % chamber pressure (if provided)
    
    % Mass
    rocket.hpcs.mass{stage, 1} = ...
        pchip(rocket.mass{stage,1}(:, 1), ...
              rocket.mass{stage,1}(:, 2));
    % Thrust
    rocket.hpcs.thrust{stage, 1} = ...
        pchip(rocket.motor.profiles.thrust.curve{stage,1}(:, 1), ...
              rocket.motor.profiles.thrust.curve{stage,1}(:, 2));
    % Burn Depth
    rocket.hpcs.burnDepth{stage,1} = ...
        pchip(rocket.motor.profiles.burnDepth.curve{stage,1}(:, 1), ...
              rocket.motor.profiles.burnDepth.curve{stage,1}(:, 2));
    % Chamber pressure
    if (flags.exists.profile.chamberPressure)
        rocket.hpcs.chamberPressure{stage,1} = ...
            pchip(rocket.motor.profiles.chamberPressure.curve{stage,1}(:, 1), ...
                  rocket.motor.profiles.chamberPressure.curve{stage,1}(:, 2));
    end
end

%% Area
% Convert given diameters (easily measurable, but not used directly in the
% process of obtaining the solution) to circular areas (not so easily
% measurable, but used directly in the process of obtaining the solution)

% Reference area for drag [m2]
rocket.area.reference = (pi/4)*rocket.cylinder.OD.^2; 

% Area of the throat for propulsion [m2]
rocket.nozzle.TA = (pi/4)*rocket.nozzle.TD.^2;

% Area of nozzle exit for adjusting thrust [m2]
rocket.nozzle.EA = (pi/4)*rocket.nozzle.ED.^2;

% Flat areas of the drogue and main parachutes when laid and fully spread
% out on the ground [m2]
rocket.parachute.drogue.area.flat = (pi/4)*rocket.parachute.drogue.flatDiameter.^2;
rocket.parachute.main.area.flat = (pi/4)*rocket.parachute.main.flatDiameter.^2;

% Inflated (hemisphere) parachute area for drag on descent [m2]
rocket.parachute.drogue.area.hemisphere = rocket.parachute.drogue.area.flat/2;
rocket.parachute.main.area.hemisphere = rocket.parachute.main.area.flat/2;

%% Local position at the launch site

% Obtain local Earth parameters
earth.launchsite.GPSh = fastinterp2(earth.terrain.longitudes, ...
                             earth.terrain.geodeticLatitudes, ...
                             earth.terrain.WGS84ToTerrain, ...
                             earth.launchsite.longitude, ...
                             earth.launchsite.latitude);
earth.launchsite.MSLh = fastinterp2(earth.terrain.longitudes, ...
                             earth.terrain.geodeticLatitudes, ...
                             earth.terrain.GeoidToTerrain, ...
                             earth.launchsite.longitude, ...
                             earth.launchsite.latitude);
[earth.time.launch.ERA, ...
 earth.time.launch.GMST, ...
 earth.time.launch.LST] = ...
    getRotAngsfromJDUT1(earth.time.launch.JD, ...
                        earth.launchsite.longitude);

% Square the rate of angular rotation
earth.pars.w2 = earth.pars.w^2;

% Obtain atmospheric density at the launch site according to the indicated 
% amosphere model
earth.time.launch.UTCYear = earth.time.launch.UTC.Year;
earth.time.launch.UTCDOY = day(earth.time.launch.UTC, 'dayofyear');
earth.time.launch.UTCSOD = earth.time.launch.UTC.Hour*3600 + ...
                           earth.time.launch.UTC.Minute*60 + ...
                           earth.time.launch.UTC.Second;
[earth.launchsite.density, ~] = atmos(earth.atmos.model, ...
                                      earth.launchsite.MSLh, ...
                                      earth.launchsite.latitude, ...
                                      earth.launchsite.longitude, ...
                                      earth.time.launch.UTCYear, ...
                                      earth.time.launch.UTCDOY, ...
                                      earth.time.launch.UTCSOD);

% Compute the Earth's rotation rate vector with respect to the ENV frame's
% orientation centered at the Earth's core
earth.vectors.w_ecef = [0; 0; earth.pars.w];
earth.T.env_ecf = ...
    getTransformationECF2ENVCoordinates(earth.launchsite.longitude, ...
                                        earth.launchsite.latitude);
earth.vectors.w_env = earth.T.env_ecf*earth.vectors.w_ecef;

% Obtain the inverse rotation for future's sake
earth.T.ecf_env = earth.T.env_ecf';

% Determine the ECEF position of the launch site based on the given inputs
earth.launchsite.PVR = computePrimeVerticalRadius(earth.pars.Req, earth.pars.e, ...
                                      earth.launchsite.latitude, "geodetic");
tmp_PVR_plus_GPSh = earth.launchsite.PVR + earth.launchsite.GPSh;
[tmp_xLaunch_ecf, tmp_yLaunch_ecf, tmp_zLaunch_ecf] = ...
    TransformGeodetic2GeocentricCoordinates(earth.launchsite.GPSh, ...
                                            earth.launchsite.latitude, ...
                                            earth.launchsite.longitude, ...
                                            earth.pars.Req, ...
                                            earth.pars.e);
earth.vectors.launchsite_ecf = [tmp_xLaunch_ecf; tmp_yLaunch_ecf; tmp_zLaunch_ecf];
clear tmp_*

%% ODE Initial conditions
% Define ODE initial conditions while rocket is resting on the launch rail
% before motor ignition. The orientation is taken with respect to the ECI
% (inertial) frame for problem solving. As such, the body-fixed angular
% rate is relative to the same ECI frame (but expressed in the body-fixed
% frame)

% First, get the necessary rotation matrices 
earth.T.rail_dcv = getTransformationR2(rocket.initial.railToVertical, "degrees");
earth.T.dcv_env = getTransformationR3(rocket.initial.eastToDownrange, "degrees");
earth.T.atLaunch.ecf_eci = getTransformationECI2ECFCoordinates(earth.time.launch.JD);

% Combine transformations
earth.T.rail_env = earth.T.rail_dcv*earth.T.dcv_env;
earth.T.atLaunch.env_eci = earth.T.env_ecf*earth.T.atLaunch.ecf_eci;

% ==================================================================
% Rotation matrix that transforms from the inertial ECI frame to the
% initial condition "rail" frame, which determines the initial conditions
% of the quaternion
earth.T.atLaunch.rail_eci = earth.T.rail_env*earth.T.atLaunch.env_eci;

% Convert this rotation matrix to an equivalent Euler-Rodrigues vector
rocket.initial.rodrigues = dcm2rod(earth.T.atLaunch.rail_eci);

% Convert this Euler-Rodrigues vector to a quaternion (the quaternion's
% scalar component is returned as the 1st element; it will be switched to
% be the 4th element while the vector component remains the same)
rocket.initial.quaternion = rod2quat(rocket.initial.rodrigues)';
% Rearrange the quaternion so that the scalar part is last
rocket.initial.quaternion = ...
    [rocket.initial.quaternion(2:4, 1); rocket.initial.quaternion(1, 1)];

% Also get the initial position of the vehicle's mass center as an
% offset away from the base of the launch rail
% 
% Prepare the rotation matrix transforming from the rail frame to the ENV
% frame (relevant since the rocket is initially aligned with the rail frame
% which makes determining the center of mass's location in the rail frame
% easy)
earth.T.env_rail = earth.T.rail_env';

% Determine the initial mass center position from the ENV origin (base of
% the launch rail) expressed in ENV coordinates
rocket.initial.position = earth.T.env_rail*[rocket.initial.railOffset; 0; ...
                   rocket.length(1,1) - rocket.wet.CoMatTip_bor{1,1}(1,1)];

% Provide a perturbation from completely completely stationary for
% numerical stability
rocket.initial.velocity = earth.T.env_rail*[0; 0; rocket.initial.speed];

% Initially, the rocket is stationary (neither translating nor rotating)
rocket.initial.angularVelocity = zeros(3, 1);

% ====================== DEFINE THE INITIAL STATE =========================
rocket.initial.state = [rocket.initial.position;
                        rocket.initial.velocity;
                        rocket.initial.quaternion;
                        rocket.initial.angularVelocity];
% =========================================================================

% Additional Rotations
earth.T.rail_env = earth.T.env_rail';
earth.T.env_dcv = earth.T.dcv_env';
earth.T.dcv_rail = earth.T.rail_dcv';
earth.T.atLaunch.eci_ecf = earth.T.atLaunch.ecf_eci';
earth.T.atLaunch.eci_env = earth.T.atLaunch.env_eci';
earth.T.atLaunch.eci_rail = earth.T.atLaunch.rail_eci';
earth.T.atLaunch.ecf_env = earth.T.atLaunch.ecf_eci*earth.T.atLaunch.eci_env;
earth.T.atLaunch.ecf_rail = earth.T.atLaunch.ecf_eci*earth.T.atLaunch.eci_rail;
earth.T.atLaunch.env_ecf = earth.T.atLaunch.ecf_env';
earth.T.atLaunch.env_rail = earth.T.atLaunch.env_ecf*earth.T.atLaunch.ecf_rail;
earth.T.atLaunch.rail_ecf = earth.T.atLaunch.ecf_rail';
earth.T.atLaunch.rail_env = earth.T.atLaunch.env_rail';


% Reorder the transformation structure T to be alphabetical
earth.T = orderfields(earth.T);
earth.T.atLaunch = orderfields(earth.T.atLaunch);

% Set options for ODE solver
flags.options.ode.settings = ...
          odeset('reltol', flags.options.ode.relativeTolerance, ...
                 'abstol', flags.options.ode.absoluteTolerance, ...
                 'initialstep', flags.options.ode.initialStep, ...
                 'maxstep', flags.options.ode.maxStep, ...
                 'stats', flags.options.ode.stats, ...
                 'refine', flags.options.ode.refine, ...
                 'events', @odevents);
clear stage