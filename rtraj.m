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
% Flags
flags.options.ode.stats = 'on'; % 'on'/'off' - Give ODE solver statistics
flags.options.ode.refine = 3; % Positive integer - Number of points to evaluate between time steps
flags.options.ode.relativeTolerance = 1e-5; % Positive scalar - Relative error tolerance bound
flags.options.ode.absoluteTolerance = 1e-7; % Positive scalar - Absolute error tolerance bound
flags.options.ode.initialStep = 1e-3; % Positive scalar - Initial step size in time
flags.options.ode.maxStep = inf; % Positive scalar - Maximum permissible time step (including inf)
flags.options.include.perturbations = false; % true/false - Consider perturbations to nominal flight
flags.options.include.cloneDynamics = false; % true/false - Save time by performing no additional calculations
flags.options.show.plots = false; % true/false - Show plots
flags.options.show.runtime = true; % true/false - Show runtime after finishing
flags.options.show.ODEtime = true; % true/false - Show solver times at each time step
flags.options.show.tablePhases = 1:6; % i:j - Sequences to display in tables (associated with DispTable)
flags.options.show.plotsPhases = 1:6; % i:j - Sequences to display in plots (associated with Plots)
flags.options.use.cache = true; % true/false - Utilize a cache to load the environment faster if applicable

%% Inputs
defineEarth % Parameters defining a 3D Earth model
defineRocket % Parameters defining the 3D flight vehicle

%% SCRIPTS 
% (Scripts see the entire workspace and add to/changes it without
% explicitly indicating the outputs)
% Earth model - gravity & magnetic fields, terrain, and atmosphere
loadEarthModel 
% Propulsion model - thrust, mass flow rate & distribution, and 
% (optional) chamber pressure
loadPropulsionData 
% Corrections and additional parameters
finalizeInputs

%% Flight Model
% Model the trajectory of a 3D multi-staged rocket launching from a ground
% facility at some designated longitude and (geodetic) latitude. The rocket
% may be suborbital (sounding) or orbital. 
% 
% Results are given with respect to the noninertial ENV frame whose origin
% is defined at the approximate GPS altitude (according to EGM2008) of the
% ground facility (above the WGS84 ellipsoid) at the launch pad's base and
% has the following axis orientations.
%   x: Points from the launch origin to true East as defined on the surface
%      of the reference ellipsoid at the given latitude and longitude of
%      the ground facility (completes the right-handed coordinate system)
%   y: Points from the launch origin to True North (NOT magnetic North) on
%      the surface of the reference ellipsoid at the given latitude and
%      longitude of the ground facility
%   z: Points from the launch origin directly upwards according to the
%      geometry of the reference ellipsoid
% 
% Orbital flights may be detected yielding further results transformed into
% the globally inertial coordinate reference frame designated J2000. Note
% that J2000 and GCRF are NOT the same coordinate systems, but are treated
% as such to avoid the transformation.
% 
% The J2000 coordinate reference frame is defined as the coordinate system
% whose origin is located at the Earth's center and has the following axis
% orientations.
%   x: Points from the Earth's center along the equator towards the vernal
%      equinox precisely as it was on January 1, 2000 12:00:00.000 TT
%   y: Points from the Earth's center along the equator into the Bay of
%      Bengal (completes the right-handed coordinate system)
%   z: Points from the Earth's center through the North Pole
% 
% The J2000 coordinates may be transformed into an Earth-fixed reference
% frame according to IAU-2000/2006 reduction using the P03 precession
% model.
[t, x, te, xe, ie] = deal(cell(rocket.stages, 1));
for s = 1:rocket.stages
    [t{s}, x{s}, te{s}, xe{s}, ie{s}] = ...
        ode113(@(tdum, xdum) odeval(tdum, xdum, earth, rocket, flags, s), ... Integrator
               [0, 2000], ... Time range
               rocket.initial.state, ...
               flags.options.ode.settings);
    t0 = t{s}(end);
    x0 = x{s}(end, :);
end