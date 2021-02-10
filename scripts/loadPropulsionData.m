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
                flags.simulate.MonteCarlo_each = false(rocket.stages, 4);
            otherwise
                rethrow(error_failedToLoadPreviousMonteCarloFlags)
        end
    end
else
    % Begin by default assuming that no profiles are to be used with a
    % Monte Carlo analysis
    flags.simulate.MonteCarlo_each = false(rocket.stages, 4);
end
%% Load thrust profile(s)

if (flags.load.profiles.thrust)
    disp("Loading thrust profile...")
    % Create temporary variable for brevity
    tmp_load = loadTheseProfiles.thrust;
    % Define thrust profile structure
    rocket.motor.profiles.thrust.name(tmp_load, 1) = "";
    rocket.motor.profiles.thrust.curve(tmp_load, 1) = {NaN};
    rocket.motor.profiles.thrust.ambientPressure(tmp_load, 1:2) = NaN;
    
    % Read each necessary file
    for stage = loadTheseProfiles.thrust
        % Define the current profile
        file = rocket.motor.files.thrust(stage, 1);
        % Read the file
        [rocket.motor.profiles.thrust.name(stage, 1), ...
         rocket.motor.profiles.thrust.curve{stage, 1}, ...
         flags.simulate.MonteCarlo_each(stage, 1), ...
         rocket.motor.profiles.thrust.ambientPressure(stage, :)] ...
            = readProfile(file);
        
        % Leave errors to Matlab regarding well-orderedness of the time
        % column
    end
else
    disp("Hitting thrust profile cache...")
end
fprintf("Loaded thrust profile\n\n")

%% Obtain mass flow rate profile(s)

if (flags.load.profiles.massFlowRate)
    disp("Loading mass flow rate profile...")
    % Create temporary variable for brevity
    tmp_load = loadTheseProfiles.massFlowRate;
    % Define structure
    rocket.motor.profiles.massFlowRate.name(tmp_load, 1) = "";
    rocket.motor.profiles.massFlowRate.curve(tmp_load, 1) = {NaN};
    
    % Read each necessary file
    for stage = loadTheseProfiles.massFlowRate
        % Define the current profile
        file = rocket.motor.files.massFlowRate(stage, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it needs to use a more simplified (linear)
        % model for this stage's behavior in mass flow rate
        if (isempty(file) || strcmp(file, ""))
            % Mark that a mass flow rate profile for this stage wasn't
            % provided, indicating that the mass of the motor will be
            % assumed to vary LINEARLY from the stage's initial mass to the
            % stage's final mass, where the final mass is simply this
            % stage's initial mass less this stage's motor mass.
            flags.exists.profile.massFlowRate(stage, 1) = false;
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [rocket.motor.profiles.massFlowRate.name(stage, 1), ...
             rocket.motor.profiles.massFlowRate.curve{stage, 1}] ...
                = readProfile(file);
            
            % Determine if a CSV profile was provided by checking the
            % values of the profile. If they're still NaN, then perform a
            % linear treatment of the mass flow rate.
            if (any(isnan(rocket.motor.profiles.massFlowRate.curve{stage, 1}), [1, 2]))
                flags.exists.profile.massFlowRate(stage, 1) = false;
            else
                % Mark that the mass flow rate profile (for this stage)
                % exists
                flags.exists.profile.massFlowRate(stage, 1) = true;
            end
        end
    end
else
    disp("Hitting mass flow rate cache...")
end
fprintf("Loaded mass flow rate profile\n\n")

%% Obtain burn depth profile(s)

if (flags.load.profiles.burnDepth)
    disp("Loading burn depth profile...")
    % Create temporary variable for brevity
    tmp_load = loadTheseProfiles.burnDepth;
    % Define structure
    rocket.motor.profiles.burnDepth.name(tmp_load, 1) = "";
    rocket.motor.profiles.burnDepth.curve(tmp_load, 1) = {NaN};
    
    % Read each necessary file
    for stage = loadTheseProfiles.burnDepth
        % Define the current profile
        file = rocket.motor.files.burnDepth(stage, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it needs to use a more simplified (linear)
        % model for this stage's behavior in mass distribution
        if (isempty(file) || strcmp(file, ""))
            % Mark that a burn-depth profile for this stage wasn't
            % provided, indicating that the mass distribution of the motor
            % will be assumed to vanish from the rocket's centerline to the
            % outside walls in radial fashion (ie BATES grain but with no
            % inner diameter - used for expression of the inertia tensor)
            flags.exists.profile.burnDepth(stage, 1) = false;
            
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [rocket.motor.profiles.burnDepth.name(stage, 1), ...
             rocket.motor.profiles.burnDepth.curve{stage, 1}] ...
                = readProfile(file);
            
            % Check if a CSV profile was defined by seeing if the func is
            % still NaN. If so, then don't do the analysis
            if (any(isnan(rocket.motor.profiles.burnDepth.curve{stage, 1}), [1, 2]))
                flags.exists.profile.burnDepth(stage, 1) = false;
            else
                % Mark that the burn depth profile (for this stage) exists
                flags.exists.profile.burnDepth(stage, 1) = true;
            end
        end
    end
else
    disp("Hitting burn depth cache...")
end
fprintf("Loaded burn depth profile\n\n")

%% Obtain chamber pressure profile(s)

if (flags.load.profiles.chamberPressure)
    disp("Loading chamber pressure profile...")
    % Create temporary variable for brevity
    tmp_load = loadTheseProfiles.chamberPressure;
    % Define structure
    rocket.motor.profiles.chamberPressure.name(tmp_load, 1) = "";
    rocket.motor.profiles.chamberPressure.curve(tmp_load, 1) = {NaN};
    
    % Read each necessary file
    for stage = loadTheseProfiles.chamberPressure
        % Define the current profile
        file = rocket.motor.files.chamberPressure(stage, 1);
        
        % Determine if the file name was empty (ie not provided) so that
        % rtraj can know if it should ignore the nozzle analysis for this
        % stage
        if (isempty(file) || strcmp(file, ""))
            % Mark that a pressure profile wasn't given for this stage,
            % indicating that any nozzle analysis should be ignored for
            % this stage
            flags.exists.profile.chamberPressure(stage, 1) = false;
        else
            % Don't define behavior if the file is spelled incorrectly such
            % that Matlab can't read it - leave those errors to Matlab
            %
            % Read the file
            [rocket.motor.profiles.chamberPressure.name(stage, 1), ...
             rocket.motor.profiles.chamberPressure.curve{stage, 1}] ...
                = readProfile(file);
            
            % Check if a CSV profile was defined by seeing if the func is
            % still NaN. If so, then don't do the analysis
            if (any(isnan(rocket.motor.profiles.chamberPressure.curve{stage, 1}), [1, 2]))
                flags.exists.profile.chamberPressure(stage, 1) = false;
            else
                % Mark that the chamber pressure profile (for this stage)
                % exists
                flags.exists.profile.chamberPressure(stage, 1) = true;
            end
        end
    end
else
    disp("Hitting chamber pressure cache...")
end
fprintf("Loaded chamber pressure profile\n\n")

%% Determine if a Monte-Carlo analysis is requested/required
flags.simulate.MonteCarlo = any(flags.simulate.MonteCarlo_each, [1, 2]);

if (flags.options.use.cache)
    motor = rocket.motor;
    % Save these results to the cache
    save(pathToCache, 'motor', 'flags', '-append')
end

% Clear up workspace of temporary variables
clear file loadTheseProfiles pathToCache previous* tmp_* stage motor