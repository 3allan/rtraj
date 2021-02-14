function f = odeval(t, x, earth, rocket, flags, stage)
% 
% Matt Werner (m.werner@vt.edu) - Dec 10, 2020
% 
% Evaluate the (highly) nonlinear function f(t, x) at time t and 
% state x for the dynamic state equation
%                        .
%                        x = f(t, x),
% 
% where x is the state determining the rocket's position, velocity, 
% orientation, and angular velocity. This process involves evaluating
% expensive terms (like gravitational potential, drag and lift
% coefficients, etc.) to get the correct expressions in evaluating the
% dynamic model of the vehicle. The state itself, x, is the 13-dimensional
% vector containing information regarding the vehicle's position, velocity,
% orientation, and angular velocity. The parameterization of the vehicle's
% orientation is taken with respect to the initial orientation of the
% launch rail, which is completely aligned with the body-fixed frame of the
% vehicle at the time of ignition and remains so until the vehicle has
% entered free flight. This orientation parameterization is realized
% through the quaternion. Explicitly, the state is defined
%         _     _
%    x = |   x   |
%        |       |
%        |   y   |
%        |       |
%        |   z   |
%        |-------|
%        |   .   |
%        |   x   |
%        |       |
%        |   .   |
%        |   y   |
%        |       |
%        |   .   |
%        |   z   |
%        |-------|
%        |  q    |
%        |   1   |
%        |       |
%        |  q    |
%        |   2   |
%        |       |
%        |  q    |
%        |   3   |
%        |       |
%        |  q    |
%        |   4   |
%        |-------|
%        |  w    |
%        |   1   |
%        |       |
%        |  w    |
%        |   2   |
%        |       |
%        |  w    |
%        |_  3  _|.
% 
% Here, the positions x(1:3) and velocities x(4:6) are taken with respect
% to the local ENV frame defined with its origin at the launch location.
% The (unit-normalized) quaternion x(7:10) tracks the body-fixed frame's
% orientation with respect to the launch rail's frame as it was at
% ignition. Finally, the angular velocity x(11:13) measures the rate of
% rotation of the body-fixed frame with respect to the same initial rail 
% frame as it was at the time of ignition. This freezing of the launch rail
% frame at the initial time makes this frame an inertial frame with respect
% to the flight vehicle. That is, this rail frame frozen at the initial
% time follows with the body, sharing the same origin as the body-fixed
% frame, serves as the frame which the body is rotating with respect to.
% 
%    Inputs:
% 
%                 t - Flight time since the initial ignition of the motor
%                     causing the rocket/flight vehicle to overcome its own
%                     weight and begin traversing up the rail guide.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%                 x - State of the rocket/flight vehicle at the current
%                     time. The state represents the current ENV
%                     position/velocity, quaternion (orientation parameter)
%                     and body-frame rotation rate expressed in the body
%                     frame.
%                     Size: 13-by-1 (vector)
%                     Units: (3-by-1) m (meters)
%                            (3-by-1) m/s (meters per second)
%                            (4-by-1) - (N/A)
%                            (3-by-1) rad/s (radians per second)
% 
%              pars - Conglomeration of input arguments related to various
%                     rocket parameters. Rocket parameters required mainly 
%                     include those related to determining the drag
%                     coefficient.
%                     Size: ?
%                     Units: ?
% 
%    Outputs:
% 
%                 f - Dynamics of the model.
%                     Size: 13-by-1 (vector)
%                     Units: (3-by-1) m/s (meters)
%                            (3-by-1) m/s2 (meters per second)
%                            (4-by-1) - (N/A)
%                            (3-by-1) rad/s2 (radians per second)
%                            
% 

% No checks

% Display integration time
if (flags.options.show.ODEtime)
    disp("t = " + t)
end

%% Setup

% Provide symbols to common elements
r_env = x(1:3, 1);
v_env = x(4:6, 1);
quat = x(7:10, 1);
w_bod = x(11:13, 1);

% Extract necessary elements
% ----------------------------------------------
% Ellipsoid equatorial radius (semimajor axis)
Req = earth.pars.Req;
% Ellipsoid eccentricity
e = earth.pars.e;
% Ellipsoid squared angular rate of rotation
w2 = earth.pars.w2;
% Angular velocity of Earth's rotation along ENV axes
w_env = earth.vectors.w_env;
% Position of launch site relative to the ECF frame
rLaunch_ecf = earth.vectors.launchsite_ecf;
% Height of launch site relative to the geoid
% Rotation matrices
Tenv_ecf = earth.T.env_ecf;
Tecf_env = earth.T.ecf_env;
% Atmosphere model
atmosModel = earth.atmos.model;
localDensity = earth.launchsite.density;
tLaunchUTC = earth.time.launch.UTC;
% Gravitational potential coefficients
nG = earth.gravity.coefficients.degrees;
mG = earth.gravity.coefficients.orders;
Cnm = earth.gravity.coefficients.cosine;
Snm = earth.gravity.coefficients.sine;

% Mass profile via HPC (during burning)
massHPC = rocket.hpcs.mass{stage,1};
%%%%%%%%%%%%%%%%%%Quarantined%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial masses of the whole stage and motor/fuel
massInitial = rocket.wet.mass(stage,1);
massMotor = rocket.motor.mass(stage,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Retarded burn time for current stage
retardedPropTime = rocket.motor.retardedTimes(stage,1);
% Burn times
burnTimes = rocket.motor.burnTime(stage,1);
% Thrust profile via HPC
thrustHPC = rocket.hpcs.thrust{stage,1};
% Exit area
areaExit = rocket.nozzle.EA(stage,1);
% Drag reference area
areaRefer = rocket.area.reference(stage,1);
% ----------------------------------------------

%% Time and Position

% ============================== Time ==================================
% Obtain current time (initial time + flight time)
tCurrentUTC = tLaunchUTC + t/86400;
% Calculate current Julian Date
tCurrentUTC_JD = juliandate(tCurrentUTC);
% Get current (D)ay (O)f (Y)ear (UTC)
tCurrentDOY = day(tCurrentUTC, 'dayofyear');
% Get current (S)econds (O)f the (D)ay since 00:00:000 UTC
tCurrentSOD = tCurrentUTC.Hour*86400 + ...
              tCurrentUTC.Minute*1440 + ...
              tCurrentUTC.Second;
          
% =========================== ECF Position ===============================
% Obtain ECF position of vehicle
r_ecf = rLaunch_ecf + Tecf_env*r_env;
% Distribute results into standard ECF (x, y, z) coordinates
xecf = r_ecf(1, 1);
yecf = r_ecf(2, 1);
zecf = r_ecf(3, 1);

% ======================== Spherical Position ============================
% Obtain spherical representation of the flight vehicle's current position
% with respect to the ECF frame
rmagecf = norm(r_ecf);
longitude_ecf = atan2(r_ecf(2), r_ecf(1));
geocentricLat_ecf = asin(r_ecf(3) / rmagecf);
geocentricColat_ecf = pi/2 - geocentricLat_ecf;

% =========================== GPS Position ===============================
% Obtain (ellipsoidal) representation of the flight vehicle's current 
% position with respect to the ECF frame associated with the earth's
% ellipsoid
[longitudeCurrent, geodeticLatCurrent, GPShCurrent] = ...
    TransformGeocentric2GeodeticCoordinates(xecf, yecf, zecf, Req, e);

% % =========================== MSL Altitude ===============================
geoidUndulationCurrent = fastinterp2(earth.terrain.longitudes, ...
                                     earth.terrain.geodeticLatitudes, ...
                                     earth.terrain.WGS84ToGeoid, ...
                                     longitudeCurrent, geodeticLatCurrent);
MSLhCurrent = GPShCurrent - geoidUndulationCurrent;

%% Atmosphere

% ========================== Atmosphere Model ============================
% Obtain atmospheric conditions at altitude
[altDensity, altTemperature] = atmos(atmosModel, MSLhCurrent, geodeticLatCurrent, longitudeCurrent, tCurrentUTC.Year, tCurrentDOY, tCurrentSOD);
% Use ideal gas law to get the pressure (unchecked, use a different atmos
% model maybe) (totally wrong - must use molecular T, not environmental T,
% and also R isn't constant)
altPressure = 287*altDensity*altTemperature;

%% Nominal Forces & Accelerations

% ========================== Transformations =============================
% Define the rotation matrix from the ECI frame to the ECF frame
Tecf_eci = getTransformationECI2ECFCoordinates(tCurrentUTC_JD);
% Define the rotation matrix from the ECF frame to the ECI frame
Teci_ecf = Tecf_eci';
% 
% Define the rotation matrix from the spherical basis set up at the ECF
% coordinates (x, y, z) back into the ECF-oriented frame
Tecf_ecfsphr = getTransformationSPHR2CARTCoordinates(longitude_ecf, geocentricColat_ecf);
% 
% Define the rotation matrix from the spherical basis set up at the ECF
% coordinates (x, y, z) to the local ENV-oriented frame
Tenv_ecfsphr = Tenv_ecf*Tecf_ecfsphr;
% 
% Define the rotation matrix from the body-fixed frame to the inertial ECI
% frame using the 4 (four) quaternions (q1, q2, q3, q4) used within the
% 13-by-1 state vector x.
Tbod_eci = getTransformationRef2Body(quat);
% 
% Define the rotation matrix from the body-fixed frame to the noninertial
% ECF frame via composition of previously defined rotation matrices
Tbod_env = Tbod_eci*Teci_ecf*Tecf_env;
Tenv_bod = Tbod_env';
% 
% Define the rotation matrix from the velocity frame to the body-fixed
% frame using the angle of attack and "velocity roll" angle
v_bod = Tbod_env*v_env; V = norm(v_bod);
[angleOfAttack, angleBod2VelRot3Roll] = calculateAngleOfAttack(v_bod, V);
Tbod_vel = getTransformationVel2Body(angleOfAttack, angleBod2VelRot3Roll);

% ==================== Gravitational Acceleration ========================
% Obtain spherical harmonic expression for the gravitational acceleration
% experienced at the current position (ECF spherical basis) (EGM2008)
accelGrav_ecfsphr = ...
    GradGravityPotential(3986004.415e8, 6378136.3, ... EGM2008 values
                         nG, mG, ... Degrees and orders
                         Cnm, Snm, ... Cosine and sine coefficients
                         rmagecf, longitude_ecf, geocentricLat_ecf);
% Rotate the expression for gravitational acceleration from the ECF
% spherical frame to the local launch ENV frame
accelGrav_env = Tenv_ecfsphr*accelGrav_ecfsphr;
% Obtain the gradient of the centrifugal potential (ECF basis)
accelCent_ecf = w2*[xecf; yecf; 0];
% Apply the exponential rate of decay in atmospheric density relative to
% the surface density at the launch site and transform the expression into
% the local launch ENV frame
densityRat = min(max(0, altDensity/localDensity), 1);
accelCentModified_env = Tenv_ecf*(densityRat*accelCent_ecf);
% State the gravitational acceleration with respect to the ENV frame
g_env = accelGrav_env + accelCentModified_env;


% ============================== Thrust ==================================
% Most of the trajectory is uncontrolled without burning any propellant;
% initially assume that the vehicle is not burning and then check if it
% actually is
flag_isFiring = false;
FT_env = zeros(3, 1);
% If not firing, then the vehicle is coasting during this stage
expelledMass = massMotor;
finalMass = massInitial - expelledMass;
mass = finalMass;
% Define the time during which the rocket is producing thrust to propel
% itself upwards
tburnCurrent = max(0, t - retardedPropTime);
% Check if the vehicle is actually firing
if (tburnCurrent < burnTimes)
    % Update flag
    flag_isFiring = true;
    % Calculate the current mass via HPC interpolation
    mass = max(mass, ppval(massHPC, tburnCurrent));
    % Reobtain the expelled mass (first used in creating the HPC - faster
    % to interpolate only once) Note: This quantity is nonpositive
    expelledMass = mass - massInitial;
    % Interpolate into the thrust profile at this current burn time
    % ensuring impossibility of numerical anomalies providing negative
    % thrust
    FT = max(0, ppval(thrustHPC, tburnCurrent) + ...
        (rocket.motor.profiles.thrust.ambientPressure(stage,1) - altPressure)*areaExit);
    % The nozzles are NOT capable of thrust vectoring (no gimbaling).
    % Therefore, the thrust is directed completely along its body-frame 3-axis
    % which points from the center of mass directly out of the nosecone, which
    % is aligned with the direction of thrust. As such, the thrust provides no
    % amount of torque on the center of mass.
    FT_bod = [0; 0; FT];
    % Transform from the body frame to the ENV frame
    FT_env = Tenv_bod*FT_bod;
end

% ================================ Drag ==================================
% Square the (noninertial) velocity of the vehicle relative to the ground
Vsq = V^2;
% Freestream dynamic pressure experienced by the vehicle at altitude
Q = 0.5*altDensity*Vsq;
% Determine the drag coefficient - Sample value (THIS IS NOT CORRECT!)
CD = 1;
% Multiply the freestream dynamic pressure by the reference area 
% >> The reference area MUST be consistent with the reference area used in
% determining the drag coefficient CD <<
Q_S = Q*areaRefer;
FD = Q_S*CD;
FD_vel = [0; 0; -FD];

% Transform from the velocity frame to the body-fixed frame
FD_bod = Tbod_vel*FD_vel;
% Transform from the velocity frame to the ENV frame
FD_env = Tenv_bod*Tbod_vel*FD_vel;

% ================================ Lift ==================================
% Sample value (THIS IS NOT CORRECT!)
FL_bod = zeros(3,1);

% Transform the lift into the env frame
FL_env = Tenv_bod*FL_bod;

%% Off-Nominal Forces & Accelerations

% ================================ Wind ==================================
% Sample value (THIS IS NOT CORRECT!)
FW_bod = zeros(3,1);

% Transform the lift into the env frame
FW_env = Tenv_bod*FW_bod;
% aoa = ?

%% Noninertial Accelerations
% Calculate contributions to inertial acceleration due to noninertial
% effects of the ENV frame (centrifugal and Coriolis forces)
accel_centrifugal = cross(w_env, cross(w_env, r_env));
accel_coriolis = 2*cross(w_env, v_env);
accel_centrifugalPlusCoriolis = accel_centrifugal + accel_coriolis;

%% Evaluate Torque Acting on the Mass Center
% ==================== Position of the Center of Mass ====================
% Center of mass is calculated relative to the nose tip in the body-fixed
% frame
CoMatTip_body = (rocket.wet.CoMatTip_body{stage,1}*rocket.wet.mass(stage,1) - ...
            rocket.motor.tipToCoM_body{stage,1}*expelledMass) / mass;


% ================== Position of the Center of Pressure ==================
% Center of pressure is calculated relative to the nose tip in the
% body-fixed frame



% ============== Aerodynamic Forces on Center of Pressure ================
FatCP_bod = FD_bod + FL_bod + FW_bod;


% ============================== Torques =================================
% Drag and lift act on the instantaneous center of pressure, creating a
% torque about the instantaneous center of mass. For rigid-body motion,
% these torques must be found in the body-fixed frame (which is very
% convenient)
% M1 = ..., M2 = ..., M3 = ...
% (Each comes very simply from expressing r x F, where r extends from the
% center of mass to the center of pressure (in the body frame) and the
% force is F = (FD_bod + FL_bod) (expressed in the body frame of course).
% That is, M1, M2, and M3 are just the body-frame components to be used in
% Euler's equation from M_bod = r_bod x F_bod
rCGtoCP_bod = [0;0;0];%[xCGtoCP_bod; yCGtoCP_bod; zCGtoCP_bod];
M_bod = cross(rCGtoCP_bod, FatCP_bod);


% ============================ Inertia Matrix =============================
Tbody_bor = getTransformationBOR2Body;
burnDepth = ppval(rocket.hpcs.burnDepth{stage,1}, retardedPropTime);
% Sample value (THIS IS NOT CORRECT!)
IGB = translateIMoIatPToCoM(mass, CoMatTip_body, ...
                            rocket.wet.IMoIatTip_body{stage,1});

% ============== Evaluate Time-Derivative of Inertia Matrix ===============
% Assume first that the rocket isn't firing any motors which is true for
% the majority of the flight. In doing so, express the 3-by-3 time
% derivative of the inertia matrix as a scalar 0 for computational
% efficiency.
dIGBdt = 0;
if (flag_isFiring)
    % Update the time derivative of the inertia matrix to reflect the fact
    % that the motor is current burning (changing shape). Thus, the mass
    % distribution is changing in time, so the derivative is nonzero.
    
    % Sample value (THIS IS NOT CORRECT!)
    dIGBdt = rocket.motor.dIMoIatCoMdt_body{stage,1};
end

%% Quaternion Kinematics
xiq = [quat(4)*eye(3) + getCrossProductEquivalentMatrix(quat(1:3)); -quat(1:3)'];

%% Nonlinear Dynamics Model
% Define the dynamics of the model - Accelerations are written in the
% noninertial ENV frame such that
% 01. f(01) = dx/dt               04. f(04) = d2x/dt2
% 02. f(02) = dy/dt               05. f(05) = d2y/dt2
% 03. f(03) = dz/dt,              06. f(06) = d2z/dt2,
% where (x, y, z) are the position components of the flight vehicle
% expressed along the ENV frame.
% 
% The kinetmatic rotation equation is 
% 07. f(07) = dq1/dt
% 08. f(08) = dq2/dt
% 09. f(09) = dq3/dt
% 10. f(10) = dq4/dt,
% where (q1, q2, q3, q4) are the quaternion components describing the
% orientation of the body-fixed frame relative to the inertial ECI frame.
% 
% Angular accelerations are written in the body-fixed frame utilizing
% torques acting about the center of mass in the body-fixed frame as
% 11. f(11) = dw1/dt
% 12. f(12) = dw2/dt
% 13. f(13) = dw3/dt,
f = zeros(13, 1);
f(1:3,   1) = v_env;
f(4:6,   1) = mass\(FT_env + FD_env + FL_env + FW_env) + g_env - accel_centrifugalPlusCoriolis;
f(7:10,  1) = 2\xiq*w_bod;
f(11:13, 1) = IGB\(M_bod - cross(w_bod, IGB*w_bod) - dIGBdt*w_bod);
