% One-time run script that calculates some additional parameters from the
% inputs of rtraj

%% Create flags
flag_ProvidedMassFlowRate = 1; % User input includes nonlinear dmdt = f(t)
flag_ProvidedChamberPressure = 1; % User input includes nonlinear pc = f(t)

%% Staging
stageNums = (1:numStages)';

%% Propulsion
% Define standard units to which the profiles must be converted
FTProfileNewUnits = ["s", "N"];
MFProfileNewUnits = ["s", "kg/s"];
PCProfileNewUnits = ["s", "Pa"];

% Correct thrust profile data from .csv file by removing excess data for
% times indicated after the burn times according to "burnTimes" and then
% ensure that the units are SI. Also ensure that the mass flow rate is
% negative.
for stage = stageNums'
    % Adjust thrust profile (necessary)
    FTProfile{stage}(FTProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
    FTProfile{stage} = convUnits(FTProfile{stage}, FTUnits__{stage}, FTProfileNewUnits);
    FTUnits__{stage} = FTProfileNewUnits;
    
    % Adjust mass flow rate profile
    if (~all(isnan(MFProfile{stage})))
        MFProfile{stage}(MFProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
        MFProfile{stage} = convUnits(MFProfile{stage}, MFUnits__{stage}, MFProfileNewUnits);
        % Ensure that the mass flow rate is negative
        MFProfile{stage}(:, 2) = -abs(MFProfile{stage}(:, 2));
    else
        % If not provided, then use linear variation
        MFProfile{stage} = -massMotor(stage)/burnTimes(stage);
        flag_ProvidedMassFlowRate = 0;
    end
    MFUnits__{stage} = MFProfileNewUnits;
    
    % Create mass profile and prepare it to be pchipped
    expelledMass = cumtrapz(MFProfile{stage}(:, 1), MFProfile{stage}(:, 2));
    % Check that the final amount doesn't exceed the specified motor mass
    excessMass = abs(expelledMass) - massMotor(stage);
    if (excessMass > 0)
        warning("Expelled mass exceeds indicated motor/fuel mass for stage $1.0f by %1.2f kg", stage, excessMass)
    end
    TMProfile{stage, 1} = [MFProfile{stage}(:, 1), massInit_(stage) + expelledMass];
    
    % Adjust chamber pressure profile
    if (~all(isnan(PCProfile{stage})))
        PCProfile{stage}(PCProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
        PCProfile{stage} = convUnits(PCProfile{stage}, PCUnits__{stage}, PCProfileNewUnits);
        PCUnits__{stage} = PCProfileNewUnits;
    else
        flag_ProvidedChamberPressure = 0;
    end
end
clear expelledMass

% Convert given diameters (easily measurable, but not used directly in the
% process of obtaining the solution) to circular areas (not so easily
% measurable, but used directly in the process of obtaining the solution)
areaRefer = pi*diaOuter_.^2/4; % Reference area for drag [m2]
areaThrt_ = pi*diaThroat.^2/4; % Area of the throat for propulsion [m2]
areaExit_ = pi*diaExit__.^2/4; % Area of nozzle exit for adjusting thrust [m2]
areaPara_ = pi*diaFlatDM.^2/8; % Inflated parachute area for drag on descent [m2]

% Obtain (H)ermite (P)olynomial (C)oefficients (HPC) for fast interpolation of
% the thrust curve, mass flow rate, and chamber pressure (if provided) per
% stage ('pchip')
% Key:
% stage, 1 - Thrust profile (Always required, so it goes first)
% stage, 2 - Mass flow rate (Optional, but still always there)
% stage, 3 - Mass (Always there)
% stage, 4 - Chamber pressure (NaN if no chamber pressure is indicated)
HPC = cell(max(stageNums), 4);
for stage = stageNums'
    % Obtain (P)iecewise (C)ubic (H)ermite (I)nterpolating (P)olynomial
    % (PCHIP) coefficients 
    %
    % PCHIP the thrust profile
    HPC{stage, 1} = pchip(FTProfile{stage}(:, 1), FTProfile{stage}(:, 2));
    %
    % PCHIP the mass flow rate
    HPC{stage, 2} = pchip(MFProfile{stage}(:, 1), MFProfile{stage}(:, 2));
    %
    % PCHIP the mass
    HPC{stage, 3} = pchip(TMProfile{stage}(:, 1), TMProfile{stage}(:, 2));
    %
    % PCHIP the chamber pressure, if provided
    if (flag_ProvidedChamberPressure)
        HPC{stage, 4} = pchip(PCProfile{stage}(:, 1), PCProfile{stage}(:, 2));
    else
        % Otherwise, give it NaN so that attempting to interpolate upon
        % HPC{..., 3} will result in a Matlab error with ppval
        HPC{stage, 4} = NaN;
    end
end
% Rename for inserting into the table 'pars'
FThpc = HPC(:, 1);
MFhpc = HPC(:, 2);
TMhpc = HPC(:, 3);
PChpc = HPC(:, 4);

% Precompute delays from launch regarding burn time
retardedPropTime = NaN(numStages, 1);
for stage = stageNums'
    if (stage == 1)
        % No delays during ignition since t = 0 occurs at ignition
        retardedPropTime(1, 1) = 0;
    else
        % Provide a nonzero offset due to burntimes, separation delays, and
        % ignition delays of previous stages
        retardedPropTime(stage, 1) = burnTimes(1, 1:stage-1) + delaySep_(1, 1:stage-1) + delayIgn_(1, 1:stage-1);
    end
end


%% Earth
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
% First, need another rotation matrix
Tenv_dcv = Tdcv_env';
% Define the distance of the mass center from the base of the rocket before
% ignition
LenNozzleCM = LenRocket(1, 1) - LenNoseCM(1, 1);

% Rearrange the quaternion so that the scalar part is last
q0 = [quat_atLaunch(2:4, 1); quat_atLaunch(1, 1)];
% Determine the initial mass center position from the ENV origin (base of
% the launch rail) expressed in ENV coordinates
r0 = Tenv_dcv*[Rail2Rckt; 0; LenNozzleCM];

% Initially, the rocket is stationary (neither translating nor rotating)
[v0, w0] = deal(zeros(3, 1));
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
Tdcv_rail = Trail_dcv';
Tenv_rail = Tenv_dcv*Tdcv_rail;
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
pars.options.dynamics = table(HavePrtrb, CloneData, ShowPlots, Runtime__, ODEtime__);
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
                    diaOuter_, diaThroat, diaExit__, diaFlatDM, ...
                    areaRefer, areaThrt_, areaExit_, areaPara_, ...
                    massInit_, massMotor, Yexhaust_, burnTimes, ...
                    delaySep_, delayIgn_, retardedPropTime, deployAlt);
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
pars.flags = table(flag_ProvidedMassFlowRate, flag_ProvidedChamberPressure);
% -----------------------------------------------------------------------