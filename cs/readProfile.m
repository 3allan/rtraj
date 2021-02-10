function [profile, func, flag_MonteCarlo, varargout] = readProfile(file)
% 
% Matt Werner (m.werner@vt.edu) - Jan 3, 2021
% 
% Extract information provided within a profile input file specifically
% formatted to be compatible with rtraj. The information contained within
% this file must include AT LEAST the vehicle's name/identification, what
% the profile to be provided actually is, and the stage during which the
% profile is valid, and the profile itself.
% 
% - File extension must be text (.txt)
% - File encoding must be UTF-8 without a BOM (byte order mark)
% 
% - '#' as the first character of the line indicates a comment
% - '>' as the first character of the line indicates an input that is NOT
%   the profile
% - '' as the first character of the line indicates the profile input (the
%   profile has no preceeding special character on every line) and its
%   entries must be separated by the comma (',')
% 
% ======================== Summary of key symbols ========================
% Comment: #
% Input: >
% Assignment: =
% Units: ()
% Pairing: ,
% Uncertainty: ?
% Relative Uncertainty: %
% 
% ================================= Rules ================================
%  - Comments:
%  1.0 # is used as the FIRST character of a line to indicate that the line
%        is to be used as a comment.
%  1.1 An empty line is treated as a comment.
% 
%  - Spacing:
%  2.0 Spaces before the FIRST character are not read (a line of spaces is
%        equivalent to an empty line).
%  2.1 Spaces on either side of an assignment (=) and commas (,) are not
%        read to arbitrary distance away.
% 
%  - Input:
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
%  3.4 , is ONLY used to provide the (csv) profile. The profile is a comma-
%        separated list with N rows and 2 columns, where the left column is 
%        populated with times (in units of seconds) and the right column is
%        some function of time evaluated at the correspondingly paired
%        times in the profile.
% 
%  - Estimation:
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
% ========================================================================
% 
%    Inputs:
% 
%              file - Path to the file containing the profile.
%                     Size: ?
%                     Units: N/A
% 
%    Outputs:
% 
%           profile - Provides a (very) short description of the intended
%                     purpose that this file has along with units for the
%                     function part (2nd column) of the provided profile,
%                     if any. Unspecified units are interpreted that the
%                     profile is unitless.
%                     Size: 1-by-1 (string)
%                     Units: ?
% 
%              func - Time-history (in seconds) of the quantity whose data
%                     is to be interpolated at various times throughout the
%                     interval. The profile is read from comma-separated
%                     value (CSV) pairs of entries, where each entry is a
%                     real number.
%                     Size: n-by-2 (matrix)
%                     Units: ?
% 
%              flag - Indicator that a Monte-Carlo simulation is requested
%                     through the use of the ? operator.
%                     Size: 1-by-1 (boolean)
%                     Units: N/A
% 
%           
% 

% USER: Skip the first line if the encoding includes a BOM (byte order
% mark). Use regular UTF-8 (without a BOM) to avoid issues in the first
% line of the file

% SCRIPT - Provides some default values in case one is requested from the
% input files (uncertainty is expected to be used)
getDefaultValues

% Open the file
fid = fopen(file, 'r', 'n', 'UTF-8');

%% Checks
% Check that the file isn't empty
thisLine = fgetl(fid);
if (thisLine == -1), error("File %s is empty", file), end
% Reset to the beginning
frewind(fid)

%% Preparation
% Set default flag to indicate that a Monte-Carlo simulation is not
% requested
flag_MonteCarlo = false;

% Get the number of lines in the file
totalLineCount = getNumberOfLines(fid);

% Begin a counter to track the current line number
ctr = 0; % Begin at 0 since the pointer is at the beginning

% Allocate memory for the profile
func = NaN(totalLineCount+2, 2);

% Initiate counter for current number of variable outputs
varargoutctr = 0;

%% Begin reading
isEndOfFile = feof(fid);
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
        
        % Ensure that the first element of VAR is not numeric
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
        
        % Handle ? wildcard by appending the default value or ABSOLUTE (no
        % %) uncertainty (whichever is requested) along with a flag
        % indicating whether a Monte-Carlo simulation is proper or not
        vals = strsplit(val, '?');
        numelVals = numel(vals);
        
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
                    error("Feature not yet implemented")
                    % Flag for Monte-Carlo is false by default
                else
                    error("Feature not yet implemented")
                    flag_MonteCarlo = true;
                end
            elseif (emptyVals(1) && ~emptyVals(2))
                % egs '> X = ?1' (absolute uncertainty)
                %     '> X = ?1%' (relative uncertainty)
                % Default value requested with specified uncertainty
                % (Monte-Carlo requested)
                requestingRelativeUncertainty = strcmp(vals{2}(end), '%');
                error("Feature not yet implemented")
                flag_MonteCarlo = true;
            elseif (~emptyVals(1) && emptyVals(2))
                % egs '> X = 1?'
                %     '> X = 1??'
                % Specified value given and requested default uncertainty
                % (Monte-Carlo requested)
                vals{2} = strcat(num2str(DEFAULT.RELATIVE_UNCERTAINTY),'%');
                VAL = computeAbsoluteUncertaintyFromChars(vals);
                flag_MonteCarlo = true;
            else
                % egs '> X = 2?1'
                %     '> X = 2?1%' 
                % Specified value given and specified uncertainty given
                % (Monte-Carlo requested)
                VAL = extractAbsoluteUncertaintyFromChars(vals);
                flag_MonteCarlo = true;
            end
        else
            % Attempt to convert the values to numeric arrays - sometimes
            % this is supposed to happen and sometimes it's not, so don't
            % take any action in the catch
            try
                VAL = cellfun(@str2num, vals);
                % Replace no uncertainty with an uncertainty of 0
                if (numelVals == 1)
                    VAL(1, 2) = 0;
                end
            catch
                % If the previous attempt to convert the values to numbers
                % failed, then the input is meant to remain as a string,
                % which, of course, may not contain a ?
                VAL = vals{1};
            end
        end
        
        % Assign specific output (vehicle, profile, stage, and func) while
        % using varargout for anything else indicated with >
        switch upper(VAR)
            case "PROFILE"
                profile = VAL;
                profileUnits = units;
            otherwise
                % Convert VAL to standard SI units
                if (isempty(units)), units = ""; end
                VAL = convUnits(VAL, units, "SI");
                % Attempt to add the input to the varargout cell. If the
                % assignment fails, then let Matlab throw the error without
                % trying to resolve it automatically
                varargout{varargoutctr + 1} = VAL;
                % Update the amount of arguments that have been assigned to
                % varargout
                varargoutctr = varargoutctr + 1;
        end
    else
        % Profile
        thisLine_timeAndFun = strsplit(thisLine, ',');
        func(ctr, :) = cellfun(@str2num, thisLine_timeAndFun);
    end
    
    % Check if this line is the end of the file
    isEndOfFile = feof(fid);
end

%% Clean up the profile
% Remove any remaining NaN values from the profile
func = func(~isnan(func));
% Reshape the profile to be n-by-2
func = reshape(func, [numel(func)/2, 2]);
% Check if there was even a profile (or if it was just a file of > inputs
% without a profile) and, if so, change its units to those in SI
if (~isempty(func))
    % Convert units (time is assumed to be in seconds)
    func(:, 2) = convUnits(func(:, 2), profileUnits, "SI");
else
    func = "No profile";
end

%% Close file
fclose(fid);