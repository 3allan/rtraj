function traj = odeden(earth, rocket, flags)
% 
% Matt Werner (m.werner@vt.edu) - Feb 13, 2021
% 
% Handle the process of numerically integrating the flight vehicle's
% dynamics prescribed by the set of ODEs for the vehicle's state. The set
% of ODEs follows the nonlinear model
%                        .
%                        x = f(t, x),
%
% where x is the state and f is the dynamics, which is dependent upon time
% t and the state. Integration begins at t = 0 when the rocket is first
% igniting and beginning to traverse the launch rail. 
% 
% At this time t = 0, the rocket is described by the initial state x(0)
% indicative of the following:
%  1. The ENV position of the vehicle's center of mass relative to the
%     launch rail's base (~[0; 0; 0])
%  2. The ENV velocity of the vehicle's center of mass relative to the
%     launch rail's base ([0; 0; 0])
%  3. The quaternion orientation parameterization of the vehicle relative
%     to the globally inertial coordinate frame (ECI) J2000
%  4. The body-fixed (principal) frame representation of the vehicle's
%    angular velocity ([0; 0; 0])
% 
% Events are created at instances of stage burnout, separation, and
% landing. Stage separation and landing warrant that integration be halted.
% 
%    Inputs:
% 
%             earth - Parameters defining an Earth model. The model
%                     includes quantities defining a reference ellipsoid,
%                     but is mainly defined by EGM2008 using high
%                     degree-order spherical harmonics for its magnetic/
%                     gravitational field models and terrain model. The
%                     atmosphere model is chosen from a list of possible
%                     options including NRLMSISE08. Transformations between
%                     several reference geodetic reference frames are
%                     included for convenience.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 
%            rocket - Parameters defining the launch vehicle. The vehicle
%                     is built from 3D components called "bodies of
%                     revolution" that are rigidly connected. The nosecone
%                     is specifiable and mass density functions are allowed
%                     to be customly defined as a function of the
%                     coordinate along the body's longitudinal axis. The
%                     motor/fuel is time-varying during the burn time.
%                     Interpolatable quantities have been pre-interpolated
%                     using PCHIP. The reference area (used for aerodynamic
%                     analysis) is the largest cross-section of the
%                     cylindrical casings in the current stage.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 
%             flags - Collection of options and indicators.
%                     Size: 1-by-1 (structure)
%                     Units: N/A
% 
%    Outputs:
% 
%              traj - Simulated trajectory of the vehicle's state through
%                     the provided dynamics.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 

% No checks

% Initial values
traj.initial.t(1,1) = 0;
traj.initial.x{1,1} = rocket.initial.state;

% Total number of separations
seps = rocket.stages + ~isinf(rocket.delays.separation(rocket.stages));

% Simulate the dynamics
for s = 1:seps
    
    % Numerically integrate the dynamics during stage s
    [traj.t{s,1}, ...
     traj.x{s,1}, ...
     traj.events.t{s,1}, ...
     traj.events.x{s,1}, ...
     traj.events.i{s,1}] = ...
         ode113(@odeval, ... Integrator and dynamics
               traj.initial.t(s,1)+[0, 100], ... Time range
               traj.initial.x{s,1}, ... Initial state
               flags.options.ode.settings, ... Settings
               earth, rocket, flags, s ... Parameters
               );
           
    % Record the state at the end of this stage
    traj.final.t(s,1) = traj.t{s,1}(end);
    traj.final.x{s,1} = traj.x{s,1}(end,:)';
    
    if (s < seps)
        % Update the initial conditions for the next stage
        traj.initial.t(s+1,1) = traj.final.t(s,1);
        traj.initial.x{s+1,1} = traj.x{s,1}(end, :)';
    end
    
    % 
end