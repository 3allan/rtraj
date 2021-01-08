% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Load the propulsion data, which requires specification of AT LEAST the
% thrust profile as a function of burn time. Additional specification of
% the mass flow rate is not required, but will be assumed to be a linear
% function of burn time varying between the stage's wet and dry masses if
% left undetermined. Specification of the motor's burn depth is also
% optional; unless explicitly given its characteristics, the motor is
% assumed to lose its mass from the rocket's centerline to its outer walls
% radially. The specification of the chamber pressure profile is also not
% required. Providing the chamber pressure profile, however, yields extra
% information about isentropic flow properties of the exhaust through the
% nozzle.
% 

% SCRIPT - Creates flags to determine which models may be used again from
% the rtraj.m workspace (if already populated by a previous simulation) and
% which require reloading
refPropulsionCache

% Only load the propulsion profiles that are either not currently stored in
% the rtraj.m workspace or have changed from the previous simulation

% Reset profile vehicle name
profiles.vehicle = "";

% Start a flag to indicate whether a Monte Carlo analysis is requested
% according to the inputs (a Monte Carlo analysis is requested if at least
% 1 input uses the ? operator with some specification of uncertainty)
if (flags.exists.previous.rtrajWorkspace)
    try
        flags.simulate.MonteCarlo_each = previous.flags.simulate.MonteCarlo_each;
    catch error_failedToLoadPreviousMonteCarloFlags
        % Let Matlab handle the errors if this assignment fails. In the case of
        % failure, the appropriate action is to clear the cache by deleting
        % any file in '<dir to rtraj.m>/cache/' related to the propulsion
        % cache. For complete assurance, run (from the rtraj.m path)
        % 'delete cache/previous_*.mat'. Also, reset the workspace by
        % running 'clearvars'. Automatic action of these commands is
        % undesired in case 
        switch error_failedToLoadPreviousMonteCarloFlags.identifier
            case 'MATLAB:nonExistentField'
                % Error results if the structure exists but the field does
                % not, which is possible if execution in this script halts
                % and quits before the flags are appended to the
                % 'previous_propulsion_profiles.mat' file in './cache/',
                % where '.' is the directory relative to rtraj.m.
                %
                % Print out the error report and give the recommended
                % procedure to repair the error (it will keep reocurring
                % as long as execution in this script quits before
                % appending to the cache).
                fprintf(2, "\nRecommended action: Clear cache and workspace.\n")
                rethrow(error_failedToLoadPreviousMonteCarloFlags)
            case 'MATLAB:undefinedVarOrClass'
                % Error results if the structure does not exist, which is
                % possible if not using the cache (obviously). Therefore,
                % define default behavior.
                flags.simulate.MonteCarlo_each = false(numStages, 4);
            otherwise
                rethrow(error_failedToLoadPreviousMonteCarloFlags)
        end
    end
else
    % Begin by default assuming that no profiles are to be used with a
    % Monte Carlo analysis
    flags.simulate.MonteCarlo_each = false(numStages, 4);
end
%% Load thrust profile(s)
if (flags.load.profiles.thrust)
    disp("Loading thrust profile...")
    % Create temporary variable for brevity
    loadtmp = loadTheseProfiles.thrust;
    % Define thrust profile structure
    profiles.thrust.vehicle(loadtmp, 1) = "";
    profiles.thrust.profile(loadtmp, 1) = "";
    profiles.thrust.stage(loadtmp, 1) = NaN;
    profiles.thrust.func(loadtmp, 1) = {NaN};
    profiles.thrust.ambientPressureAtTest(loadtmp, [1,2]) = NaN;
    profiles.thrust.motorMass(loadtmp, [1,2]) = NaN;
    profiles.thrust.Yexhaust(loadtmp, [1,2]) = NaN;
    
    % Read each necessary file
    for ii = loadTheseProfiles.thrust
        % Define the current profile
        file = profiles.thrust.path(ii, 1);
        % Read the file
        [ ...
            profiles.thrust.vehicle(ii, 1), ...
            profiles.thrust.stage(ii, 1), ...
            profiles.thrust.profile(ii, 1), ...
            profiles.thrust.func{ii, 1}, ...
            flags.simulate.MonteCarlo_each(ii, 1), ...
            profiles.thrust.motorMass(ii, :), ...
            profiles.thrust.Yexhaust(ii, :), ...
            profiles.thrust.ambientPressureAtTest(ii, :) ...
            ] ...
            = readProfile(file);
        % Perform some basic checks on this stage's thrust profile
        checkProfile(profiles.thrust, "THRUST", ii, file)
    end
else
    disp("Hitting thrust profile cache...")
end
fprintf("Loaded thrust profile\n\n")

% Check thrust profile for consistency between vehicle and profile names
% and if all equal, reassign the vehicle and profile fields to be a single
% string
checkProfile(profiles.thrust)
profiles.thrust.vehicle = profiles.thrust.vehicle(1);
profiles.thrust.profile = profiles.thrust.profile(1);

%% Obtain mass flow rate profile(s)
if (flags.load.profiles.massFlowRate)
    disp("Loading mass flow rate profile...")
    % Create temporary variable for brevity
    loadtmp = loadTheseProfiles.massFlowRate;
    % Define structure
    profiles.massFlowRate.vehicle(loadtmp, 1) = "";
    profiles.massFlowRate.profile(loadtmp, 1) = "";
    profiles.massFlowRate.stage(loadtmp, 1) = NaN;
    profiles.massFlowRate.func(loadtmp, 1) = {NaN};
    
    % Read each necessary file
    for ii = loadTheseProfiles.massFlowRate
        % Define the current profile
        file = profiles.massFlowRate.path(ii, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it needs to use a more simplified (linear)
        % model for this stage's behavior in mass flow rate
        if (isempty(file) || strcmp(file, ""))
            % Mark that a mass flow rate profile for this stage wasn't
            % provided, indicating that the mass of the motor will be
            % assumed to vary LINEARLY from the stage's initial mass to the
            % stage's final mass, where the final mass is simply this
            % stage's initial mass less this stage's motor mass.
            flags.exists.profile.massFlowRate(ii, 1) = false;
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [ ...
                profiles.massFlowRate.vehicle(ii, 1), ...
                profiles.massFlowRate.stage(ii, 1), ...
                profiles.massFlowRate.profile(ii, 1), ...
                profiles.massFlowRate.func{ii, 1}, ...
                flags.simulate.MonteCarlo_each(ii, 2) ...
                ] ...
                = readProfile(file);
            % Perform some basic checks on this stage's mass flow rate profile
            checkProfile(profiles.massFlowRate, "MASS-FLOW-RATE", ii, file)
            
            % Check if a CSV profile was defined by seeing if the func is
            % still NaN. If so, then don't do the analysis
            if (any(isnan(profiles.massFlowRate.func{ii, 1}), [1, 2]))
                flags.exists.profile.massFlowRate(ii, 1) = false;
            else
                % Mark that the mass flow rate profile (for this stage)
                % exists
                flags.exists.profile.massFlowRate(ii, 1) = true;
            end
        end
    end
else
    disp("Hitting mass flow rate cache...")
end
fprintf("Loaded mass flow rate profile\n\n")

% Check mass flow rate profile for consistency between vehicle and profile
% names and, if all equal, reassign the vehicle and profile fields to be a
% single string
checkProfile(profiles.massFlowRate)
profiles.massFlowRate.vehicle = profiles.massFlowRate.vehicle(1);
profiles.massFlowRate.profile = profiles.massFlowRate.profile(1);

%% Obtain burn depth profile(s)
if (flags.load.profiles.burnDepth)
    disp("Loading burn depth profile...")
    % Create temporary variable for brevity
    loadtmp = loadTheseProfiles.burnDepth;
    % Define structure
    profiles.burnDepth.vehicle(loadtmp, 1) = "";
    profiles.burnDepth.profile(loadtmp, 1) = "";
    profiles.burnDepth.stage(loadtmp, 1) = NaN;
    profiles.burnDepth.func(loadtmp, 1) = {NaN};
    
    % Read each necessary file
    for ii = loadTheseProfiles.burnDepth
        % Define the current profile
        file = profiles.burnDepth.path(ii, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it needs to use a more simplified (linear)
        % model for this stage's behavior in mass distribution
        if (isempty(file) || strcmp(file, ""))
            % Mark that a burn-depth profile for this stage wasn't
            % provided, indicating that the mass distribution of the motor
            % will be assumed to vanish from the rocket's centerline to the
            % outside walls in radial fashion (ie BATES grain but with no
            % inner diameter - used for expression of the inertia tensor)
            flags.exists.profile.burnDepth(ii, 1) = false;
            
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [ ...
                profiles.burnDepth.vehicle(ii, 1), ...
                profiles.burnDepth.stage(ii, 1), ...
                profiles.burnDepth.profile(ii, 1), ...
                profiles.burnDepth.func{ii, 1}, ...
                flags.simulate.MonteCarlo_each(ii, 3) ...
                ] ...
                = readProfile(file);
            % Perform some basic checks on this stage's mass flow rate profile
            checkProfile(profiles.burnDepth, "BURN-DEPTH", ii, file)
            
            % Check if a CSV profile was defined by seeing if the func is
            % still NaN. If so, then don't do the analysis
            if (any(isnan(profiles.burnDepth.func{ii, 1}), [1, 2]))
                flags.exists.profile.burnDepth(ii, 1) = false;
            else
                % Mark that the burn depth profile (for this stage) exists
                flags.exists.profile.burnDepth(ii, 1) = true;
            end
        end
    end
else
    disp("Hitting burn depth cache...")
end
fprintf("Loaded burn depth profile\n\n")

% Check burn depth profile for consistency between vehicle and profile
% names and, if all equal, reassign the vehicle and profile fields to be a
% single string
checkProfile(profiles.burnDepth)
profiles.burnDepth.vehicle = profiles.burnDepth.vehicle(1);
profiles.burnDepth.profile = profiles.burnDepth.profile(1);

%% Obtain chamber pressure profile(s)
if (flags.load.profiles.chamberPressure)
    disp("Loading chamber pressure profile...")
    % Create temporary variable for brevity
    loadtmp = loadTheseProfiles.chamberPressure;
    % Define structure
    profiles.chamberPressure.vehicle(loadtmp, 1) = "";
    profiles.chamberPressure.profile(loadtmp, 1) = "";
    profiles.chamberPressure.stage(loadtmp, 1) = NaN;
    profiles.chamberPressure.func(loadtmp, 1) = {NaN};
    
    % Read each necessary file
    for ii = loadTheseProfiles.chamberPressure
        % Define the current profile
        file = profiles.chamberPressure.path(ii, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it should ignore the nozzle analysis for this
        % stage
        if (isempty(file) || strcmp(file, ""))
            % Mark that a pressure profile wasn't given for this stage,
            % indicating that any nozzle analysis should be ignored for
            % this stage
            flags.exists.profile.chamberPressure(ii, 1) = false;
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [ ...
                profiles.chamberPressure.vehicle(ii, 1), ...
                profiles.chamberPressure.stage(ii, 1), ...
                profiles.chamberPressure.profile(ii, 1), ...
                profiles.chamberPressure.func{ii, 1}, ...
                flags.simulate.MonteCarlo_each(ii, 4) ...
                ] ...
                = readProfile(file);
            % Perform some basic checks on this stage's mass flow rate profile
            checkProfile(profiles.chamberPressure, "CHAMBER-PRESSURE", ii, file)
            
            % Check if a CSV profile was defined by seeing if the func is
            % still NaN. If so, then don't do the analysis
            if (any(isnan(profiles.chamberPressure.func{ii, 1}), [1, 2]))
                flags.exists.profile.chamberPressure(ii, 1) = false;
            else
                % Mark that the chamber pressure profile (for this stage)
                % exists
                flags.exists.profile.chamberPressure(ii, 1) = true;
            end
        end
    end
else
    disp("Hitting chamber pressure cache...")
end
fprintf("Loaded chamber pressure profile\n\n")

% Check chamber pressure profile for consistency between vehicle and
% profile names and, if all equal, reassign the vehicle and profile fields
% to be a single string
checkProfile(profiles.chamberPressure)
profiles.chamberPressure.vehicle = profiles.chamberPressure.vehicle(1);
profiles.chamberPressure.profile = profiles.chamberPressure.profile(1);

%% Compare names of vehicles and profiles for each profile to ensure they're
% all the same and, if so, combine them in each profile to reduce repetition
checkStringEquivalence([profiles.thrust.vehicle, ...
    profiles.massFlowRate.vehicle, profiles.burnDepth.vehicle, ...
    profiles.chamberPressure.vehicle], false)
% List the vehicle's name in its profiles from the thrust profile (the
% thrust profile MUST be present to define any rocket)
profiles.vehicle = profiles.thrust.vehicle;
% Get number of fields (to reorganize - more profiles get added later)
numelfields = numel(fieldnames(profiles));
% Reorder the fields to move the vehicle's name to the top
if (numelfields == 5)
    % Mass profile has not yet been added
    profiles = orderfields(profiles, {'vehicle', ...
        'thrust', 'massFlowRate', 'burnDepth', 'chamberPressure'});
elseif (numelfields == 6)
    % Mass profile has been added
    profiles = orderfields(profiles, {'vehicle', ...
        'thrust', 'massFlowRate', 'burnDepth', 'chamberPressure', 'mass'});
end

%% Determine if a Monte-Carlo analysis is requested/required
flags.simulate.MonteCarlo = any(flags.simulate.MonteCarlo_each, [1, 2]);

if (flags.options.use.cache)
    % Save these results to the cache
    save(pathToCache, 'profiles', 'flags', '-append')
end

% Clear up workspace of temporary variables
clear file pathToCache loadtmp previous