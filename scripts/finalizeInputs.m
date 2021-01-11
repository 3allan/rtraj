% One-time run script that calculates some additional parameters from the
% inputs of rtraj

%% Staging
stageNums = (1:rocket.stages)';

%% Propulsion

% Extract the burn times from the thrust profiles
for stage = stageNums'
    % Explicitly place the burn times of the motors in with the thrust
    % profile, where the first column indicates the nominal burn time and
    % the second column indicates the (absolute) uncertainty
    profiles.thrust.burnTimes(stage, :) = [profiles.thrust.func{stage}(end, 1), 0];
end

% Units for propulsive characteristics are already SI and burn times are
% assumed to be cut off when thrust reaches 0.
% 
% Objectives:
%  1. Provide a (linear) mass flow rate if a mass flow rate profile wasn't
%     explicitly given
%  1.i. Ensure that the mass flow rate is negative
%  2. Provide a (linear variation of) mass distribution if a burn depth
%     profile wasn't explicitly given
%  3. Indicate whether a chamber pressure profile was provided or not
for stage = stageNums'
    %% Mass Flow Rate Profile
    % Ensure that the mass flow rate profile is negative
    if (flags.exists.profile.massFlowRate(stage))
        % Negate the provided mass flow rate profile
        profiles.massFlowRate.func{stage}(:, 2) = ...
            -abs(profiles.massFlowRate.func{stage}(:, 2));
    else
        % Temporary variables for brevity and clarity
        tmp_tb = profiles.thrust.burnTimes(stage, 1);
        tmp_m0 = massInit(stage);
        tmp_massFlowRate = -tmp_m0/tmp_tb;
        % Define the linear profile
        profiles.massFlowRate.func{stage, :} = [  0,     tmp_massFlowRate
                                                tmp_tb,  tmp_massFlowRate]; 
    end
    
    %% Mass Profile
    % Temporary variables for brevity and clarity
    tmp_t = profiles.massFlowRate.func{stage}(:, 1);
    tmp_f = profiles.massFlowRate.func{stage}(:, 2);
    tmp_motorMass = profiles.thrust.motorMass(stage, 1);
    % Create mass profile from mass flow rate profile and prepare it to be pchipped
    tmp_expelledMass = -cumtrapz(tmp_t, tmp_f);
    % Check that the final amount doesn't exceed the specified motor mass
    tmp_excessMass = tmp_expelledMass(end) - tmp_motorMass;
    if (tmp_excessMass > 0.001)
        warning("Expelled mass exceeds indicated motor/fuel mass for stage %1.0f by %1.2f kg", stage, 0.001+tmp_excessMass)
    end
    % Define the mass profile
    tmp_mass = massInit(stage) + tmp_expelledMass;
    profiles.mass.func{stage, 1} = [tmp_t, tmp_mass];
    profiles.mass.expelled{stage, 1} = tmp_expelledMass;
    
    %% Burn Depth Profile
    % Nothing to do if the profile is provided, but if it's not provided,
    % then make a linear variation varying from the rocket's centerline to
    % its outer diameter0
    if (~flags.exists.profile.burnDepth)
        % Temporary variable for brevity and clarity
        tmp_ID = diaOuter(stage); % Note: ID (inner diameter) used with OD (outer diameter)
        profiles.burnDepth.func{stage, :} = [  0,         0 
                                             tmp_tb,    tmp_ID];
    end
    
    %% Chamber Pressure Profile
    % Nothing to do since the flags regarding existence of a chamber
    % pressure profile have already been determined
    
    %% Clean up temporary variables
    clear clear tmp_* expelledMass
end

%% Area
% Convert given diameters (easily measurable, but not used directly in the
% process of obtaining the solution) to circular areas (not so easily
% measurable, but used directly in the process of obtaining the solution)
area.reference = pi*diaOuter.^2/4; % Reference area for drag [m2]
area.throat = pi*diaThroat.^2/4; % Area of the throat for propulsion [m2]
area.exit = pi*diaExit.^2/4; % Area of nozzle exit for adjusting thrust [m2]
area.parachute = pi*diaFlatDM.^2/8; % Inflated parachute area for drag on descent [m2]

%% Interpolating Coefficients
% Obtain (H)ermite (P)olynomial (C)oefficients (HPC) for fast interpolation
% ('pchip') of the thrust, mass flow rate, burn depth, and chamber pressure
% (if provided) per each stage.
for stage = stageNums'
    % Temporary variables for brevity and clarity
    tmp_tthrust = profiles.thrust.func{stage}(:, 1);
    tmp_fthrust = profiles.thrust.func{stage}(:, 2);
    %
    tmp_tmassFlowRate = profiles.massFlowRate.func{stage}(:, 1);
    tmp_fmassFlowRate = profiles.massFlowRate.func{stage}(:, 2);
    %
    tmp_tmass = profiles.mass.func{stage}(:, 1);
    tmp_fmass = profiles.mass.func{stage}(:, 2);
    %
    tmp_tburnDepth = profiles.burnDepth.func{stage}(:, 1);
    tmp_fburnDepth = profiles.burnDepth.func{stage}(:, 2);
        
    % Obtain (P)iecewise (C)ubic (H)ermite (I)nterpolating (P)olynomial
    % (PCHIP) coefficients
    %
    % PCHIP the thrust profile
    HPC.thrust{stage, 1} = pchip(tmp_tthrust, tmp_fthrust);
    %
    % PCHIP the mass flow rate
    HPC.massFlowRate{stage, 1} = pchip(tmp_tmassFlowRate, tmp_fmassFlowRate);
    %
    % PCHIP the mass
    HPC.mass{stage, 1} = pchip(tmp_tmass, tmp_fmass);
    %
    % PCHIP the chamber pressure, if possible
    if (flags.exists.profile.chamberPressure)
        HPC.chamberPressure{stage, 1} = pchip(tmp_tburnDepth, tmp_fburnDepth);
    end
    
    % Clean up temporary variables
    clear tmp_*
end

%% Precompute Delays From Launch Regarding Burn Time
retardedPropTime = NaN(numStages, 2);
for stage = stageNums'
    if (stage == 1)
        % No delays during ignition since t = 0 occurs at ignition
        retardedPropTime(1, :) = zeros(1, 2);
    else
        % Provide a nonzero offset due to burntimes, separation delays, and
        % ignition delays of previous stages
        retardedPropTime(stage, :) = profiles.thrust.burnTimes(1:stage-1, :) + ...
            delaySep(1:stage-1, :) + delayIgn(1:stage-1, :);
    end
end

%% Earth

% Obtain local Earth parameters
GPShLaunchsite = fastinterp2(longitudes, geodeticLatitudes, WGS84ToTerrain, LonLaunch, LatLaunch);
MSLhLaunchsite = fastinterp2(longitudes, geodeticLatitudes, GeoidToTerrain, LonLaunch, LatLaunch);
[launchsite.ERA.initial, launchsite.GMST.initial, launchsite.LST.initial] = ...
    getRotAngsfromJDUT1(JDLaunch, deg2rad(LonLaunch));



% Square the rate of angular rotation
w2 = w^2;


% Convert indicates launch longitude and latitude to radians
LonLaunchrad = deg2rad(LonLaunch);
LatLaunchrad = deg2rad(LatLaunch);


% Obtain atmospheric density at the launch site according to the indicated 
% amosphere model
tLaunchUTCYear = tLaunchUTC.Year;
tLaunchUTCDOY = day(tLaunchUTC, 'dayofyear');
tLaunchUTCSOD = tLaunchUTC.Hour*86400 + ...
              tLaunchUTC.Minute*1440 + ...
              tLaunchUTC.Second;
[localDensity, ~] = atmos(atmosModel, MSLhLaunchsite, LatLaunchrad, LonLaunchrad, tLaunchUTCYear, tLaunchUTCDOY, tLaunchUTCSOD);


% Compute the Earth's rotation rate vector with respect to the ENV frame's
% orientation centered at the Earth's core
w_ecef = [0; 0; w];
Tenv_ecf = getTransformationECF2ENVCoordinates(LonLaunchrad, LatLaunchrad);
w_env = Tenv_ecf*w_ecef;
% Rename by transposing the above results to pass along with pars
[w_ecef_rowvec, w_env_rowvec] = deal(w_ecef', w_env');
% Obtain the inverse rotation for future's sake
Tecf_env = Tenv_ecf';


% Determine the ECEF position of the launch site based on the given inputs
RNLaunch = computePrimeVerticalRadius(Req, e, LatLaunchrad, "geodetic");
RNLaunch_plus_hLaunch = RNLaunch + GPShLaunchsite;
[xLaunch_ecf, yLaunch_ecf, zLaunch_ecf] = TransformGeodetic2GeocentricCoordinates(GPShLaunchsite, LatLaunchrad, LonLaunchrad, Req, e);
rLaunch_ecef = [xLaunch_ecf; yLaunch_ecf; zLaunch_ecf];
rLaunch_ecef_rowvec = rLaunch_ecef';
clear RNplush xLaunch_ecf yLaunch_ecf zLaunch_ecf

% %% Structures
% % Determine the location of the fully wetted vehicle's center of mass
% % relative to the body-fixed frame measured from the nose tip
% wetMassCenterPosition_Nose = nan(3, 1, numStages);
% for stage = stageNums'
%     wetMassCenterPosition_Nose = [0, 0, -LenNoseCM(stage)];
%     dryMassCenterPosition_Nose = 
% end


%% ODE Initial conditions
% Define ODE initial conditions while rocket is resting on the launch rail
% before motor ignition. The orientation is taken with respect to the ECI
% (inertial) frame for problem solving. As such, the body-fixed angular
% rate is relative to the same ECI frame (but expressed in the body-fixed
% frame)
% 
% First, get the necessary rotation matrices 
Tecf_eci_atLaunch = getTransformationECI2ECFCoordinates(JDLaunch);
Tenv_eci_atLaunch = Tenv_ecf*Tecf_eci_atLaunch;
Tdcv_env = getTransformationR3(East2DwnR, "degrees");
Trail_dcv = getTransformationR2(Rail2Vert, "degrees");
Trail_eci_atLaunch = Trail_dcv*Tdcv_env*Tenv_eci_atLaunch;
% Convert this rotation matrix to an equivalent Euler-Rodrigues vector
rod_atLaunch = dcm2rod(Trail_eci_atLaunch);
% Convert this Euler-Rodrigues vector to a quaternion (scalar part is
% returned as the first element)
quat_atLaunch = rod2quat(rod_atLaunch)';
% 
% Also get the initial position of the vehicle's mass center as an
% offset away from the base of the launch rail
% 
% First, need some more rotation matrices
Tenv_dcv = Tdcv_env';
Tdcv_rail = Trail_dcv';
Tenv_rail = Tenv_dcv*Tdcv_rail;
% Define the distance of the mass center from the base of the rocket before
% ignition
LenNozzleCM = LenRocket(1, 1) - LenNoseCM(1, 1);

% Rearrange the quaternion so that the scalar part is last
q0 = [quat_atLaunch(2:4, 1); quat_atLaunch(1, 1)];
% Determine the initial mass center position from the ENV origin (base of
% the launch rail) expressed in ENV coordinates
r0 = Tenv_rail*[Rail2Rckt; 0; LenNozzleCM];

% Initially, the rocket is stationary (neither translating nor rotating)
w0 = zeros(3, 1);
% Provide a perturbation from completely completely stationary for
% numerical stability
v0 = Tenv_rail*[0; 0; veps];
x0 = [r0; v0; q0; w0];
% Define the initial time at which integration begins
t0 = 0;
% Define the final time at which integration ends
tf = 2000;

% Set options for ODE solver
odeOpts = odeset('reltol', odereltol, 'abstol', odeabstol, ...
                 'initialstep', odeinstep, 'maxstep', odemxstep, ...
                 'stats', odestats_, 'refine', oderefine, ...
                 'events', @odevents);
             
%% Additional Rotations
Trail_env = Tenv_rail';
Teci_ecf_atLaunch = Tecf_eci_atLaunch';

%% Categorization
% -----------------------------------------------------------------------
% Place necessary variables into pars categorized by application and amount
% of rows (tables require equal amounts of rows per structure) (structures
% are accessed by . notation (Ex: pars.time.tLocalNow.TimeZone indicates
% that 'pars' has a structure called 'time' which holds a value called
% 'tLocalNow' which has a property called 'TimeZone')).
% --- ODE ---
pars.options.dynamics = table(HavePrtrb, CloneData, ShowPlots, Runtime, ODEtime);
pars.options.numericalSolvers = table(odeOpts);
% --- TIME ---
pars.time = table(tLocalNow, JDLaunch, tLaunch, tLaunchUTC, tLaunchUTCYear, tLaunchUTCDOY, tLaunchUTCSOD);
% --- LAUNCH ---
% Pass quantities to do with the launch time and initial orientation of the
% rail with respect to the ENV frame.
pars.launchsite = table(LonLaunch, LatLaunch, LonLaunchrad, LatLaunchrad, ... 
                        ERA0, GMST0, LST0, GPShLaunchsite, MSLhLaunchsite, ...
                        localTemp, localDensity, towerSpan, rLaunch_ecef_rowvec);
pars.launchrail = table(Rail2Rckt, Rail2Vert, East2DwnR, veps);
% --- FLIGHT VEHICLE ---
% (geometric/phase-defining) together
pars.rocket = table(stageNums, LenRocket, LenNoseCM, ...
                    diaOuter , diaThroat, diaExit  , diaFlatDM, ...
                    areaRefer, areaThrt , areaExit , areaPara , ...
                    massInit , massMotor, Yexhaust , burnTimes, ...
                    delaySep , delayIgn , retardedPropTime, deployAlt);
% --- PROPULSION ---
pars.thrustProfile = table(FTProfile, FThpc); % thrust
pars.massFlowRate = table(MFProfile, MFhpc); % mass flow rate
pars.mass = table(TMProfile, TMhpc); % mass
pars.chamberPressure = table(PCProfile, PChpc); % chamber pressure
% --- AERODYNAMICS ---
pars.dragCoefficient = table(RasAeroCD); % drag coefficient
% --- EARTH MODEL ---
pars.ellipsoid = table(GM, Req, Rpo, f, e, w_ecef_rowvec, w_env_rowvec, w, w2); % ellipsoid parameters
pars.terrain.longitudes = table(longitudes); % longitudes of height maps
pars.terrain.geodeticLatitudes = table(geodeticLatitudes); % geodetic latitudes of height maps
pars.terrain.heights = table(WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain); % height maps
pars.coefficients.gravity = table(nG, mG, Cnm, Snm); % gravitational field potential harmonic coefficients
pars.coefficients.magnetf = table(nM, mM, gnm, hnm); % magnetic field potential harmonic coefficients
pars.atmosphericModel = table(atmosModel); % atmosphere model
% --- Rotation Matrices ---
pars.rotations.back.singles = table(Trail_dcv, Tdcv_env, Tenv_ecf, Tecf_eci_atLaunch);
pars.rotations.forw.singles = table(Teci_ecf_atLaunch, Tecf_env, Tenv_dcv, Tdcv_rail);
pars.rotations.forw.jumps = table(Tenv_rail);
pars.rotations.back.jumps = table(Trail_env);
% --- FLAGS ---
pars.flags = table(flags);
% -----------------------------------------------------------------------