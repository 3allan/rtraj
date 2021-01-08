function checkProfile(profile, profileName, stage, file)
% 
% Matt Werner (m.werner@vt.edu) - Jan 6, 2021
% 
% Ensure that particular fields of the provided profile are consistent with
% the expected behavior. Expected behavior of a profile, which is a
% structure object containing information about the particular profile (eg
% thrust, mass flow rate, etc.), is determined from the respective and
% provided input files.
% 
% The expected behaviors are:
%   1) All vehicle names within the "vehicle" field must be identical to
%      one another.
%   2) All profile names within the "profile" field must be identical to
%      one another AND match the profile name (with spaces " " filled by
%      hyphens "-"). 
%      Profile names are:
%       - Thrust profile: "THRUST"
%       - Mass Flow Rate profile: "MASS-FLOW-RATE"
%       - Burn Depth profile: "BURN-DEPTH"
%       - Chamber Pressure profile: "CHAMBER-PRESSURE"
%   3) CSV pairs of time t and profile function f(t), resulting in a time-
%      history matrix with many rows (t, f(t)), exhibit strictly monotonic
%      increasing behavior in time t. (e.g., time = [0; 0.1; 2; 30.5; 34; ...])
% 
%    Inputs:
% 
%           profile - Provided profile whose fields are to be checked for
%                     consistency.
%                     Size: 1-by-1 (struct)
%                     Units: N/A
% 
%       profileName - Optional(!) Name of the profile ("THRUST", etc.).
%                     Size: 1-by-1 (string)
%                     Units: N/A
% 
%             stage - Optional(!) Current stage number in which to check
%                     the profile. If left unprovided, then the profile is
%                     checked for consistency in its fields "vehicle" and
%                     "profile". Otherwise, the profile is checked that the
%                     read stage number from the file matches the indicated
%                     stage, that the "profile" field matches the intended
%                     profile, and that time is strictly monotonically
%                     increasing.
%                     Size: 1-by-1 (scalar)
%                     Units: N/A
% 
%              file - Optional(!) Path to the file from which the profile
%                     was read (used to point to source if a check fails,
%                     resulting in an error).
%                     Size: 1-by-1 (string)
%                     Units: N/A
%           

% Enforce that there are only 1 or 4 inputs
if (1 < nargin && nargin < 4)
    error("Invalid number of inputs (must be 1 or 4)")
end

% Determine behavior of the check based on number of inputs
if (nargin == 1)
    % Check that each string field within the fields of the structure are
    % exactly equivalent
    numelProfileVehicles = numel(profile.vehicle);
    if (numelProfileVehicles == 1)
        % No need to compare values to themselves
        return
    end
    % Compare elements in fields to each other
    strVehicle = profile.vehicle(1, 1);
    strProfile = profile.profile(1, 1);
    for ii = 2:numelProfileVehicles
        % Check for mismatching of vehicle names
        if (~strcmp(strVehicle, profile.vehicle(ii, 1)))
            error("Vehicle mismatch.")
        end
        % Check for mismatching of profile purposes
        if (~strcmp(strProfile, profile.profile(ii, 1)))
            error("Profile mismatch.")
        end
    end
    return
end
        
% Basic checks for data within file
if (profile.stage(stage, 1) ~= stage)
    issueReadError("Profile does not match stage", file, -1)
end
if (~strcmp(profile.profile(stage, 1), profileName))
    issueReadError("Profile does not match profile", file, -1)
end
if (~all(diff(profile.func{stage, 1}(:, 1)) > 0))
    issueReadError("Profile time must be monotonic increasing", file, -1)
end