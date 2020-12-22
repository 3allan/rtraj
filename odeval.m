function f = odeval(t, x, pars)
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

%% Setup

% Provide symbols to common elements
r_env = x(1:3, 1);
v_env = x(4:6, 1);
quat = x(7:10, 1);
w_body = x(11:13, 1);

x_L = r_env(1); % Longitudinal ENV coordinate (L)ambda (east)
x_P = r_env(2); % Geodetic latitudinal ENV coordinate (P)hi (north)
x_H = r_env(3); % Ellipsoidal height ENV coordinate (H) (up)

% Extract necessary elements
% ----------------------------------------------
% Ellipsoid equatorial radius (semimajor axis)
Req = pars.ellipsoid.Req;
% Ellipsoid eccentricity
e = pars.ellipsoid.e;
% Angular velocity of Earth's rotation along ENV axes
w_env = pars.launchsite.w_env_rowvec';
% Position of launch site relative to the ECF frame
rLaunch_ecf = pars.launchsite.rLaunch_ecef_rowvec';
% Rotation matrices
Tenv_ecf = pars.rotations.Tenv_ecf;
Tecf_env = pars.rotations.Tecf_env;
% Atmosphere model
atmosModel = pars.atmosphericModel.atmosModel;
tLaunchUTC = pars.time.tLaunchUTC;
% Gravitational potential coefficients
nG = pars.coefficients.gravity.nG;
mG = pars.coefficients.gravity.mG;
Cnm = pars.coefficients.gravity.Cnm;
Snm = pars.coefficients.gravity.Snm;
% Thrust profile via HPC
FThpc = pars.thrustProfile.FThpc;
% ----------------------------------------------

%% Intermediate calculations

% Obtain current time (initial time + flight time)
tCurrentUTC = tLaunchUTC + t/86400;
% Get current (D)ay (O)f (Y)ear (UTC)
tCurrentDOY = day(tCurrentUTC, 'dayofyear');
% Get current (S)econds (O)f the (D)ay since 00:00:000 UTC
tCurrentSOD = tCurrentUTC.Hour*86400 + tCurrentUTC.Minute*1440 + tCurrentUTC.Second;

% Obtain ECF position of vehicle
r_ecf = rLaunch_ecf + Tecf_env*r_env;
x = r_ecf(1, 1);
y = r_ecf(2, 1);
z = r_ecf(3, 1);

% Obtain ellipsoidal representation of the flight vehicle's current
% position with respect to the ECF frame
[longitudeCurrent, geodeticLatCurrent, GPShCurrent] = TransformGeocentric2GeodeticCoordinates(x, y, z, Req, e);
% Obtain atmospheric conditions
[d, T] = atmos(atmosModel, GPShCurrent, geodeticLatCurrent, longitudeCurrent, tCurrentUTC.Year, tCurrentDOY, tCurrentSOD);
% Use ideal gas law to get the pressure (unchecked, use a different atmos model maybe)
p = 287*d*T;

%% Nominal Forces
% Compute the nominal forces acting on the flight vehicle (weight, lift,
% thrust, drag). Each force is rather involved in its expression, depending
% generally upon time, but also other effects like body-geometry, the
% atmosphere, angle of attack, etc.

% Obtain spherical representation of the flight vehicle's current position
% with respect to the ECF frame
rmagecf = norm(r_ecf);
lonecf = atan2(r_ecf(2), r_ecf(1));
latecf = asin(r_ecf(3) / rmagecf);
colatecf = pi/2 - latecf;

% Gravity (function doesn't yet exist)
FG_geocentricSpherical = GradGravityPotential(nG, mG, Cnm, Snm, rmagecf, lonecf, latecf, GM, Req);
% Obtain transformation matrices
Tcart_sphr = getTransformationSPHR2CARTCoordinates(lonecf, colatecf);
% Transform from ECF spherical coordinates to ENV using a rotation
FG_env = Tenv_ecf*Tecf_sph*FG_geocentricSpherical;

% Thrust
FT = ppval(FThpc, t);
FT_bod = [0, 0, FT]';
% Obtain transformation matrices
% Trail_bod = % quaternion rotation matrix
% Tenv_rail = % easy rotation (in document)
% Transform from the body frame to the ENV frame
FT_env = Tenv_bod*FT_bod;

% Drag
% FD = ????;
FD_vel = [0, 0, -FD];
% Transform from the velocity frame to the ENV frame
FD_env = Tenv_bod*Tbod_vel*FD_vel;

% Lift

%% Noninertial accelerations
% Calculate contributions to inertial acceleration due to noninertial
% effects of the ENV frame (centrifugal and Coriolis forces)
accel_centrifugal = cross(w_env, cross(w_env, r_env));
accel_coriolis = 2*cross(w_env, v_env);
accel_centrifugalPlusCoriolis = accel_centrifugal + accel_coriolis;

%% Nonlinear Dynamics Model
% Define the dynamics of the model - Accelerations are written in the
% noninertial ENV frame such that
% 1. f(1) = dx/dt                4. f(4) = d2x/dt2
% 2. f(2) = dy/dt                5. f(5) = d2y/dt2
% 3. f(3) = dz/dt                6. f(6) = d2z/dt2,
% where (x, y, z) are the position components of the flight vehicle
% expressed along the ENV frame.
f = zeros(13, 1);
f(1:3, 1) = v_env;
f(4:6, 1) = ... 
                - accel_centrifugalPlusCoriolis;
