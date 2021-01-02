% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Cache system to improve the load time of performance-based profiles (eg
% thrust profile) of the flight vehicle.
% 

% Set default load behavior
loadTheseFTProfiles = 1:numStages;
loadTheseMFProfiles = 1:numStages;
loadThesePCProfiles = 1:numStages;

% Set default flags indicating that no cache currently exists and the
% propulsion profiles determining motor/engine characteristics should be
% loaded as if for the first time
% 
% Make use of flag_rtraj_workspace_populated from Earth cache reference and
% assume that the propulsion profiles exist as well
flags.exists.previous.PropulsionCache = false;
flags.load.ThrustProfile = true;
flags.load.MassFlowRateProfile = true;
flags.load.ChamberPressureProfile = true;

% Document the last modified-by date for the profiles
FTModDate = getModifiedByDate(FTprofile); % Thrust
MFModDate = getModifiedByDate(MFprofile); % Mass flow rate
PCModDate = getModifiedByDate(PCprofile); % Chamber pressure

% Attempt to load previous file name for the propulsion model into the workspace.
% This process is justified because the load only occurs once per run of
% rtraj.m and the resulting table contains few entries, so the load is fast
try
    % Obtain the path that points to where the propulsion cache is
    pathToCache = getPathToCache('previous_propulsion_profiles.mat');
    % Attempt to load the cache (may not exist)
    previous_propulsion_profiles = load(pathToCache);
    flags.exists.previous.PropulsionCache = true;
catch error_CacheMiss
    % Check possible causes for error
    switch error_CacheMiss.identifier
        case 'MATLAB:load:couldNotReadFile'
            % Create a new cache file since one currently doesn't exist and
            % carry on to load the propulsion profiles
            save(pathToCache, 'FTprofile', 'FTModDate', ...
                'MFprofile', 'MFModDate', 'PCprofile', 'PCModDate')
        otherwise
            rethrow(error_CacheMiss)
    end
    return
end

% Already known if the workspace from a previous simulation still exists or
% not. Assume that individual elements of the workspace were not deleted so
% that the variables may be used again if the inputs didn't change



if (flags.exists.previous.rtrajWorkspace)
    % Determine which propulsion profiles need to be loaded/updated and which
    % may be reused from the previous simulation
    %
    % Check which, if any, thrust profiles differ from the previous
    % simulation
    previous_FTprofile = previous_propulsion_profiles.FTprofile;
    previous_FTModDate = previous_propulsion_profiles.FTModDate;
    loadTheseFTProfiles = hasFileChanged(FTprofile, FTModDate, ...
                                previous_FTprofile, previous_FTModDate);
    clear previous_FTprofile previous_FTModDate
    %
    % Check which, if any, mass flow rate profiles differ from the previous
    % simulation
    previous_MFprofile = previous_propulsion_profiles.MFprofile;
    previous_MFModDate = previous_propulsion_profiles.MFModDate;
    loadTheseMFProfiles = hasFileChanged(MFprofile, MFModDate, ...
                                previous_MFprofile, previous_MFModDate);
    clear previous_MFprofile previous_MFModDate
    %
    % Check which, if any, chamber pressure profiles differ from the previous
    % simulation
    previous_PCprofile = previous_propulsion_profiles.PCprofile;
    previous_PCModDate = previous_propulsion_profiles.PCModDate;
    loadThesePCProfiles = hasFileChanged(PCprofile, PCModDate, ...
                                previous_PCprofile, previous_PCModDate);
end

% Adjust load flags according to any changes made in the files
if (isempty(loadTheseFTProfiles)), flags.load.ThrustProfile = false; end
if (isempty(loadTheseMFProfiles)), flags.load.MassFlowRateProfile = false; end
if (isempty(loadThesePCProfiles)), flags.load.ChamberPressureProfile = false; end

% Update the cache
save(pathToCache, 'FTprofile', 'FTModDate', ...
    'MFprofile', 'MFModDate', 'PCprofile', 'PCModDate')