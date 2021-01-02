function [vehicle, profile, stage, mass, ID, OD, Y, p, FT, flag] = ...
    readThrustProfile(file)
% 
% Matt Werner (m.werner@vt.edu) - Jan 1, 2021
% 
% Extract information provided within a propulsion file specifically
% formatted to be compatible with rtraj. The information contained within
% this file pertains ONLY the motor and not the surrounding body of the
% rocket in any way. As such, mass refers to the (initial) motor mass
% and, likewise, the inner-diameter (ID) and outer-diameter (OD) refer to
% the (mean) inner diameter and the outer diameter of the solid motor as
% measured from the motor's centerline, respectively. The (mean) ID is
% requested to be used in determining the inertia matrix of the solid motor
% while burning, which is assumed to be a BATES grain (hollow cylinder with
% a cross-section consisting of two (2) concentric circles) for the
% purpose of rotational dynamics.
% 
% ======================== Summary of key symbols ========================
% Comment: #
% Input: >
% Assignment: =
% Units: ()
% Pairing: ,
% Estimate: ?
% Percent: % (4.0d only)
% 
% ================================= Rules ================================
%  Comments:
%  1.0 # is used as the FIRST character of a line to indicate that the line
%        is to be used as a comment.
%  1.1 An empty line is treated as a comment.
% 
%  Spacing:
%  2.0 Spaces before the FIRST character are not read (a line of spaces is
%        equivalent to an empty line).
%  2.1 Spaces on either side of an assignment (=) and commas (,) are not
%        read to arbitrary distance away.
% 
%  Input:
%  3.0 > is used as the FIRST character of a line to indicate that the line
%        is to be used for assigning a parameter. 
%  3.1 Such a line beginning with > expects an assignment (=) to exist.
%  3.2 Any parameter whose first character is a number is assumed to be a
%        number and NOT a string. Thus, input meant to be read as text may
%        NOT begin with a number.
%  3.3 (...) is used as the LAST set of characters in a line for parameter
%        input to indicate that parameter's units, where ... represents the
%        units (like 'm', 'kg', 's', 'N', etc. Also recognizes and accepts 
%        (some) standard units (like 'ft', 'lb', 's', 'lbf', etc.) Note the
%        distinction between lb (pounds, weight) and lbf (pounds, force).
%  (Thrust profile)
%  3.4 , is ONLY used to provide the thrust profile. The thrust profile is
%        a comma-separated list with N rows and 2 columns, where the left
%        column is populated with burn times (in units of seconds) of the
%        motor and the right column indicates the associated thrust of the
%        motor at some ambient pressure.
% 
% Estimation
%  4.0 ? is used SUCCEEDING AN ASSIGNMENT (=) to indicate that the input is
%        not precisely known and an estimate is to be used. 
%        a) When used alone (eg X = ?), this estimate is constant and
%           does NOT give reason to do repeated simulations. Effectively,
%           the assignment X = ? is
%                               X = DEFAULT_VALUE,
%           where the units are taken to be SI.
%        b) When used with an immediately preceeding number (eg X = 1?),
%           a default amount of uncertainty is automatically applied and a
%           Monte-Carlo type analysis begins performing a specified amount
%           of iterations. Effectively, the assignment X = Y? is 
%                     X = Y +/- DEFAULT_RELATIVE_UNCERTAINTY,
%           where Y is a known, specified number.
%        c) When used with immediately preceeding and succeeding numbers 
%           (eg X = 1?2), then the succeeding number is treated as the 
%           ABSOLUTE uncertainty and a Monte-Carlo type analysis begins
%           performing a specified amount of iterations. Effectively, 
%           the assignment X = Y?Z is
%                                 X = Y +/- Z,
%           where Y and Z are known, specified numbers of the same units.
%        d) When used with immediately preceeding and succeding numbers
%           followed by % (eg X = 1?2%), then the succeeding number is
%           treated as the RELATIVE uncertainty and a Monte-Carlo type
%           analysis begins performing a specified amount of iterations.
%           Effectively, the assignment X = Y?Z% is
%                                X = Y (1 +/- Z),
%           where Y and Z are known, specified numbers.
%        Note that units, provided by (), may still succeed all of the
%        above cases.
% 
%    Inputs:
% 
%              file - Path to the thrust profile file that provides
%                     information regarding motor characteristics.
%                     Size: ?
%                     Units: N/A
% 
%    Outputs: vehicle, profile, stage, mass, ID, OD, Y, FT
% 
%           vehicle - Name of the rocket.
%                     Size: ?
%                     Units: N/A
% 
%           profile - Indicates the type AND units of profile being
%                     supplied/prepared to read. For this particular
%                     file, the profile should be assigned "THRUST"
%                     followed by the corresponding units that the thrust
%                     values carry.
%                     Size: 1-by-1 (string)
%                     Units: ?
% 
%             stage - Stage number during which this motor is actively
%                     burning.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%              mass - Initial mass of the motor before ignition.
%                     Size: 1-by-1 (scalar)
%                     Units: kg (kilograms)
% 
%                ID - Inner diameter of the motor as measured from its
%                     centerline.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                OD - Outer diameter of the motor as measured from its
%                     centerline.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 Y - Exhaust ratio of specific heats.
%                     Size: 1-by-1 (scalar)
%                     Units: - (unitless)
% 
%              flag - Indicator that a Monte-Carlo simulation is requested
%                     through the use of the ? operator
% 

% Set default flag to indicate that a Monte-Carlo simulation is not
% requested
flag = false;

% Open the file
fid = fopen(file, 'r', 'n', 'UTF-8');

% Check if file is empty
isEndOfFile = feof(fid);
if (isEndOfFile)
    error("Thrust profile is empty.")
end

% Get the number of lines in the file
totalLineCount = getNumberOfLines(fid);

% Reset the head back to the beginning since it was just moved to the end
frewind(fid);

% Skip the first line to avoid issues with the byte order mark (BOM) and
% check if the first line was the file's only line
thisLine = fgetl(fid);
isEndOfFile = feof(fid);
if (isEndOfFile)
    error("No thrust profile provided.")
end

% Begin a counter to track the current line number
ctr = 1; % Begin at 1 since the first line was skipped

% Allocate memory for the thrust profile
FT = NaN(totalLineCount, 2);

% Begin reading
while (~isEndOfFile)
    % Get the current line and increment line counter
    thisLine = fgetl(fid);
    ctr = ctr + 1;
    
    % Remove all spaces (fast)
    thisLine = strrep(thisLine, ' ', '');
    
    % Skip the line if it is blank (or was filled with some spaces) or if
    % the first character indicates that the line is a comment
    if (isempty(thisLine) || strcmp('#', thisLine(1)))
        isEndOfFile = feof(fid);
        continue
    end
    
    % Process the line
    if (strcmp('>', thisLine(1)))
        % Current line is a parameter input with an assignment
        % Split current line at the assignment operator (=)
        assignment=thisLine(2:end);
        splitline = strsplit(assignment, '=', 'CollapseDelimiters', false);
        % This cell could be:
        % - a 1x1 cell with an empty char (e: line is blank after '>')
        % - a 1x1 cell with the original string (e: line is missing '=')
        % - a 1x2 cell with the original string split where '=' is
        %   - 1st element is an empty char (e: no variable)
        %   - 2nd element is an empty char (e: no value)
        %   - Both elements are populated (goal)
        % - a 1xn cell with the original string split where '=' is, but
        %   more than once (n times) (e: multiple =)
        numelsplitline = numel(splitline);
        if (numelsplitline == 1)
            % No assignment (eg '> X 1')
            issueReadError("No assignment to variable", file, ctr)
        elseif (numelsplitline > 2)
            % Too many assignments (eg '> X == 1' or '> X = 1 = 1', etc.)
            issueReadError("Too many assignments issued", file, ctr)
        elseif (isempty(splitline{1}))
            % Missing variable (eg '> = 1')
            issueReadError("No variable for assignment", file, ctr)
        elseif (isempty(splitline{2}))
            % Missing value (eg '> X = ')
            issueReadError("No value for assignment", file, ctr)
        end
        
        % Safe to assign the variable and value
        VAR = splitline{1};
        val = splitline{2};
        
        % Ensure that the first element of var is not numeric
        if (isstrprop(VAR(1), 'digit'))
            issueReadError("Variable cannot begin with digit", file, ctr)
        end
        
        % Use regular expression (regexp) to extract units from the value
        [valCells, unitsCell] = regexp(val, '(?<=\()(.*?)(?=\s*\))', 'split', 'match');
        % Get value (independent of any units)
        if (isempty(unitsCell))
            % No split parenthesis - whole value went into one cell
            val = valCells{1};
        else
            % Parentheses are split apart - assign the value, but exclude
            % the opening parenthesis '(', from the first cell (safe since
            % errors on input have already been checked)
            val = valCells{1}(1:end-1);
        end
        %
        numelunitsCell = numel(unitsCell);
        % Determine if the value is unitless or not, or if an error should
        % be thrown
        if (numelunitsCell == 0)
            % Unitless
            units = '';
        elseif (numelunitsCell == 1)
            % Take units out from cell
            units = convertCharsToStrings(unitsCell{1});
        else
            % Too many parenthesis pairs () issued
            issueReadError("Invalid units", file, ctr)
        end
        % Other indicators of unitless quantities (using parentheses) are:
        % '-', 'N/A', 'NA', 'NONE'
        if (~strcmp(units, '') && (strcmp(units, '-') || ...
                strcmp(units, 'N/A') ||  strcmp(units, 'NA') || ...
                strcmp(units, 'NONE')))
            units = '';
        end
        
        % Handle ? wildcard by appending the default value or ABSOLUTE
        % uncertainty (whichever is requested) along with a flag indicating
        % whether a Monte-Carlo simulation is proper or not
        vals = strsplit(val, '?');
        numelVals = numel(vals);
        % Replace an uncertainty of 0 with no uncertainty
        if (numelVals > 1 && (strcmp(vals{2}, '0') || strcmp(vals{2}, '0%')))
            % Remove ? entirely and treat the input as perfect
            vals(:, 2) = [];
        end
        % Find which entries are empty
        emptyVals = cellfun(@isempty, vals);
        numelVals = numel(vals);
        % Determine if the variable should be assigned a default value or
        % be associated with some uncertainty
        if (numelVals > 1)
            % ? operator requested - emptyVals is 1x2 logical array
            if (all(emptyVals))
                % egs '> X = ?'
                %     '> X = ??'
                % Default value requested with/without the default
                % uncertainty. In the case that the uncertainty is not
                % requested (X = ?), then no Monte-Carlo request is made.
                % Otherwise (X = ??), the Monte-Carlo request is made
                numQuestionMarks = numel(strfind(val, '?'));
                if (numQuestionMarks == 1)
                    error("Feature not yet implemented (missing default file)")
                    % Flag for Monte-Carlo is false by default
                else
                    error("Feature not yet implemented (missing default file)")
                    flag = true;
                end
            elseif (emptyVals(1) && ~emptyVals(2))
                % eg '> X = ?1'
                % Default value requested with specified RELATIVE uncertainty
                % (Monte-Carlo requested)
                error("Feature not yet implemented (missing default file)")
                flag = true;
            elseif (~emptyVals(1) && emptyVals(2))
                % eg '> X = 1?'
                % Specified value given and requested default uncertainty
                % (Monte-Carlo requested)
                error("Feature not yet implemented (missing default file)")
                flag = true;
            else
                % egs '> X = 2?1'
                %     '> X = 2?1%' 
                % Specified value given and specified uncertainty given
                % (Monte-Carlo requested)
                
                % Check if % is specified at the end of the input
                hasPercentOnEnd = strcmp(vals{2}(end), '%');
                if (hasPercentOnEnd)
                    % Override relative uncertainty with absolute
                    % uncertainty by first removing %, converting to
                    % arrays, and then performing the conversion to
                    % absolute uncertainty
                    vals{2} = vals{2}(1:end-1); % Remove %
                end
                VAL = cellfun(@str2num, vals); % Convert to numbers
                if (hasPercentOnEnd)
                    VAL(2) = VAL(1)*VAL(2)/100; % Convert to abs.
                end
                flag = true;
            end
        else
            % Attempt to convert the values to numeric arrays - sometimes
            % this is supposed to happen and sometimes it's not, so don't
            % take any action in the catch
            try
                VAL = cellfun(@str2num, vals);
            catch
                % If the previous attempt to convert the values to numbers
                % failed, then the input is meant to remain as a string,
                % which, of course, may not contain a ?
                VAL = vals{1};
            end
        end
        
        
        % Assign quantities based on whether the value is
        % unitless/nondimensional or dimensional
        if (isempty(units))
            % Match var name with select list of options for unitless
            % inputs
            switch VAR
                case "VEHICLE"
                    vehicle = VAL;
                case "STAGE"
                    stage = VAL;
                case "EXHAUSTHEATCAPACITYRATIO"
                    Y = VAL;
                otherwise
                    % Issue a warning that this input won't be used in
                    % simulation
                    warning("Unassigned input (%s: line %1.0f)", file, ctr)
            end
        else
            % Match var name with select list of options for dimensional
            % inputs
            switch VAR
                case "PROFILE"
                    profile = [VAL, units];
                case "MASS"
                    mass = convUnits(VAL, units, "kg");
                case "INNERDIAMETER"
                    ID = convUnits(VAL, units, "m");
                case "OUTERDIAMETER"
                    OD = convUnits(VAL, units, "m");
                case "AMBIENTPRESSURE"
                    p = convUnits(VAL, units, "Pa");
                otherwise
                    % Issue a warning that this input won't be used in
                    % simulation
                    warning("Unassigned input (%s: line %1.0f)", file, ctr)
            end
        end
    else
        % Thrust profile
        tFT_thisLine = strsplit(thisLine, ',');
        FT(ctr, :) = cellfun(@str2num, tFT_thisLine);
    end
    
    
    % Check if this line is the end of the file
    isEndOfFile = feof(fid);
end

% Remove any remaining NaN values from the thrust profile
FT = FT(~isnan(FT));
% Reshape
FT = reshape(FT, [numel(FT)/2, 2]);
% Convert units of thrust (time is assumed to be in seconds)
FT(:, 2) = convUnits(FT(:, 2), profile(1, 2), "N");
% Remove units from the profile output since it has just been converted to
% Newtons (SI)
profile(:, 2) = [];


% Close the file
fclose(fid);