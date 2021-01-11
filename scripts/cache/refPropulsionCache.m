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
    profiles.thrust.path(:, 2) = ...
        getModifiedByDate(profiles.thrust.path(:, 1));
    % Mass flow rate
    profiles.massFlowRate.path(:, 2) = ...
        getModifiedByDate(profiles.massFlowRate.path(:, 1)); 
    % Burn Depth
    profiles.burnDepth.path(:, 2) = ...
        getModifiedByDate(profiles.burnDepth.path(:, 1));
    % Chamber pressure
    profiles.chamberPressure.path(:, 2) = ...
        getModifiedByDate(profiles.chamberPressure.path(:, 1));
    
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
    pathToCache = getPathToCache('previous_propulsion_profiles.mat');
    % Attempt to load the cache (may not exist)
    previous = load(pathToCache);
    flags.exists.previous.PropulsionCache = true;
catch error_cacheMiss
    % Check possible causes for error
    switch error_cacheMiss.identifier
        case 'MATLAB:load:couldNotReadFile'
            % Create a new cache file since one currently doesn't exist and
            % carry on to load the propulsion profiles
            save(pathToCache, 'profiles')
        otherwise
            rethrow(error_cacheMiss)
    end
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
    FTprofile = profiles.thrust.path(:, 1);
    FTModDate = profiles.thrust.path(:, 2);
    previous_FTprofile = previous.profiles.thrust.path(:, 1);
    previous_FTModDate = previous.profiles.thrust.path(:, 2);
    loadTheseProfiles.thrust = hasFileChanged(FTprofile, FTModDate, ...
                                previous_FTprofile, previous_FTModDate)';
    clear FTprofile FTModDate previous_FTprofile previous_FTModDate
    %
    %
    % Check which, if any, mass flow rate profiles differ from the previous
    % simulation
    MFprofile = profiles.massFlowRate.path(:, 1);
    MFModDate = profiles.massFlowRate.path(:, 2);
    previous_MFprofile = previous.profiles.massFlowRate.path(:, 1);
    previous_MFModDate = previous.profiles.massFlowRate.path(:, 2);
    loadTheseProfiles.massFlowRate = hasFileChanged(MFprofile, MFModDate, ...
                                previous_MFprofile, previous_MFModDate)';
    clear MFprofile MFModDate previous_MFprofile previous_MFModDate
    %
    %
    % Check which, if any, burn depth profiles differ from the previous
    % simulation
    BDprofile = profiles.burnDepth.path(:, 1);
    BDModDate = profiles.burnDepth.path(:, 2);
    previous_BDprofile = previous.profiles.burnDepth.path(:, 1);
    previous_BDModDate = previous.profiles.burnDepth.path(:, 2);
    loadTheseProfiles.burnDepth = hasFileChanged(BDprofile, BDModDate, ...
                                previous_BDprofile , previous_BDModDate)';
    clear BDprofile BDModDate previous_BDprofile previous_BDModDate
    %
    %
    % Check which, if any, chamber pressure profiles differ from the previous
    % simulation
    PCprofile = profiles.chamberPressure.path(:, 1);
    PCModDate = profiles.chamberPressure.path(:, 2);
    previous_PCprofile = previous.profiles.chamberPressure.path(:, 1);
    previous_PCModDate = previous.profiles.chamberPressure.path(:, 2);
    loadTheseProfiles.chamberPressure = hasFileChanged(PCprofile, PCModDate, ...
                                previous_PCprofile, previous_PCModDate)';
    clear PCprofile PCModDate previous_PCprofile previous_PCModDate
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
save(pathToCache, 'profiles')