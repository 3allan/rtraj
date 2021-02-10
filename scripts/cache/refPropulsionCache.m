% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Cache system to improve the load time of performance-based profiles (eg
% thrust profile) of the flight vehicle.
% 

% Set default load behavior
loadTheseProfiles.thrust = 1:rocket.stages;
loadTheseProfiles.massFlowRate = 1:rocket.stages;
loadTheseProfiles.burnDepth = 1:rocket.stages;
loadTheseProfiles.chamberPressure = 1:rocket.stages;

% Set default flags indicating that no cache currently exists and the
% propulsion profiles determining motor/engine characteristics should be
% loaded as if for the first time. Also set a flag indicating that a Monte
% Carlo analysis should not be conducted
% 
% Make use of flag_rtraj_workspace_populated from Earth cache reference and
% assume that the propulsion profiles exist as well
flags.exists.previous.PropulsionCache = false;
flags.load.profiles.thrust = true;
flags.load.profiles.massFlowRate = true;
flags.load.profiles.chamberPressure = true;
flags.load.profiles.burnDepth= true;

% Check if the cache system should be used. If not, stop here; otherwise,
% continue
if (~flags.options.use.cache), return, end

% Document the last modified-by date for the profiles (each field on the
% RHS is 2-by-1 before assignment)
try
    % Thrust
    rocket.motor.files.thrust(:, 2) = ...
        getModifiedByDate(rocket.motor.files.thrust(:, 1));
    % Mass flow rate
    rocket.motor.files.massFlowRate(:, 2) = ...
        getModifiedByDate(rocket.motor.files.massFlowRate(:, 1)); 
    % Burn Depth
    rocket.motor.files.burnDepth(:, 2) = ...
        getModifiedByDate(rocket.motor.files.burnDepth(:, 1));
    % Chamber pressure
    rocket.motor.files.chamberPressure(:, 2) = ...
        getModifiedByDate(rocket.motor.files.chamberPressure(:, 1));
    
catch error_gettingModifiedByDates
    switch error_gettingModifiedByDates.identifier
        case 'MATLAB:subsassigndimmismatch'
            fprintf(2, "\nStaged data must align columnwise.\n")
            rethrow(error_gettingModifiedByDates)
        otherwise
            rethrow(error_gettingModifiedByDates)
    end
end

% Attempt to load previous file name for the propulsion model into the workspace.
% This process is justified because the load only occurs once per run of
% rtraj.m and the resulting table contains few entries, so the load is fast
try
    % Obtain the path that points to where the propulsion cache is
    pathToCache = getPathToCache('previous_propulsion_motor.mat');
    % Attempt to load the cache (may not exist)
    previous = load(pathToCache);
    flags.exists.previous.PropulsionCache = true;
catch flag_cacheMiss
    % Check possible causes for error
    switch flag_cacheMiss.identifier
        case 'MATLAB:load:couldNotReadFile'
            % Create a new cache file since one currently doesn't exist and
            % carry on to load the propulsion profiles
            motor = rocket.motor;
            save(pathToCache, 'motor')
            clear motor
        otherwise
            rethrow(flag_cacheMiss)
    end
    % Mark the flag
    flags.exists.previous.PropulsionCache = false;
    return
end

% Already known if the workspace from a previous simulation still exists or
% not from the Earth cache. Assume that individual elements of the workspace 
% were not deleted so that the variables may be used again if the inputs
% didn't change
if (flags.exists.previous.rtrajWorkspace)
    % Determine which propulsion profiles need to be loaded/updated and which
    % may be reused from the previous simulation
    %
    % Check which, if any, thrust profiles differ from the previous
    % simulation
    loadTheseProfiles.thrust = ...
        hasFileChanged(rocket.motor.files.thrust(:, 1), ...
                       rocket.motor.files.thrust(:, 2), ...
                       previous.motor.files.thrust(:, 1), ...
                       previous.motor.files.thrust(:, 2));
    %
    % Check which, if any, mass flow rate profiles differ from the previous
    % simulation
    loadTheseProfiles.massFlowRate = ...
        hasFileChanged(rocket.motor.files.massFlowRate(:, 1), ...
                       rocket.motor.files.massFlowRate(:, 2), ...
                       previous.motor.files.massFlowRate(:, 1), ...
                       previous.motor.files.massFlowRate(:, 2));
    %
    % Check which, if any, burn depth profiles differ from the previous
    % simulation
    loadTheseProfiles.burnDepth = ...
        hasFileChanged(rocket.motor.files.burnDepth(:, 1), ...
                       rocket.motor.files.burnDepth(:, 2), ...
                       previous.motor.files.burnDepth(:, 1), ...
                       previous.motor.files.burnDepth(:, 2));
    %
    % Check which, if any, chamber pressure profiles differ from the previous
    % simulation
    loadTheseProfiles.chamberPressure = ...
        hasFileChanged(rocket.motor.files.chamberPressure(:, 1), ...
                       rocket.motor.files.chamberPressure(:, 2), ...
                       previous.motor.files.chamberPressure(:, 1), ...
                       previous.motor.files.chamberPressure(:, 2));
end

% Adjust load flags according to any changes made in the files
if (isempty(loadTheseProfiles.thrust))
    flags.load.profiles.thrust = false;
end
if (isempty(loadTheseProfiles.massFlowRate))
    flags.load.profiles.massFlowRate = false;
end
if (isempty(loadTheseProfiles.burnDepth))
    flags.load.profiles.burnDepth = false; 
end
if (isempty(loadTheseProfiles.chamberPressure))
    flags.load.profiles.chamberPressure = false;
end

% Update the cache
motor = rocket.motor;
save(pathToCache, 'motor')
clear motor