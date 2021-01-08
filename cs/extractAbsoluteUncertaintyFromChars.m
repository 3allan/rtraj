function VAL = computeAbsoluteUncertaintyFromChars(vals)
% 
% Matt Werner (m.werner@vt.edu) - Jan 4, 2021
% 
% Calculate the absolute uncertainty in a value given the value and its
% relative or absolute uncertainty. Given the absolute uncertainty, no
% calculation is necessary, though conversion from character to scalar is
% required.
% 
%    Inputs:
% 
%              vals - Contains the nominal value and the relative
%                     uncertainty, where both quantities are characters.
%                     Size: 1-by-2 (cell)
%                     Units: N/A
% 
%    Outputs:
% 
%              VALS - Contains the nominal value and the absolute
%                     uncertainty, where both quantities are scalars.
%                     Size: 1-by-2 (cell)
%                     Units: N/A
% 


% Check if % is specified at the end of the input
requestingRelativeUncertainty = strcmp(vals{2}(end), '%');
if (requestingRelativeUncertainty)
    % Override relative uncertainty with absolute
    % uncertainty by first removing %, converting to
    % arrays, and then performing the conversion to
    % absolute uncertainty
    vals{2} = vals{2}(1:end-1); % Remove %
end
VAL = cellfun(@str2num, vals); % Convert to numbers
if (requestingRelativeUncertainty)
    VAL(2) = VAL(1)*VAL(2)/100; % Convert to abs.
end