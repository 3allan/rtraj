function y = convertUnits(x, units, newUnits)
% 
% Matt Werner (m.werner@vt.edu) - Dec 9, 2020
% 
% Convert the quantities contained in x, which hold the indicated units, to
% the the indicated desired units. The unit conversions are assumed to
% happen column-wise (i.e. in the "1" direction) for standard 2D matrices x.
% 
%    Inputs:
% 
%                 x - Quantity holding values of particular units that are
%                     desired to be changed. The units of x are changed in
%                     each column, i.e. column-wise.
%                     Size: n-by-m (scalar)
%                     Units: ?
% 
%             units - Specified units of each column in x.
%                     Size: 1-by-m (string)
%                     Units: - (N/A)
% 
%          newUnits - Specified units of each column in x that the provided
%                     values in x shall be converted to.
%                     Size: 1-by-m (string)
%                     Units: - (N/A)
% 

% Check that inputs are strings
checkInput(units, newUnits)

% Leave if inputs are valid to Matlab (output unassigned)

% Check that sizes are compatible
sizex = size(x);
if (~all(sizex == size(units)) || ~all(sizex == size(newUnits)))
    error("Incompatabile sizes of quantities and units.")
end

% Allocate space to hold unit conversions
unitmult = nan(1, size(x, 2));

% Determine which unit conversions go where in `unitmult'
for ii = 1:sizex
    % Lowercase inputs
    currentUnit = lower(units(ii));
    desiredUnit = lower(newUnits(ii));
    
    %% Time
%     switch currentUnit
%         case {"millisecond", "milliseconds", "ms"}
%             switch desiredUnit
    
    %% Distance
    switch currentUnit
        case {"millimeter", "millimeters", "mm"}
            switch desiredUnit
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 1;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 0.1;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 0.001;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 0.000001;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 0.0394;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 0.0033;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 6.2137e-7;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 5.3996e-7;
            end
        case {"centimeter", "centimeters", "cm"}
            switch desiredUnit
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 10;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 1;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 0.01;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 0.00001;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 0.394;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 0.033;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 6.2137e-6;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 5.3996e-6;
            end
        case {"meter", "meters", "m"}
            switch desiredUnit
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 1000;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 100;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 1;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 0.001;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 39.3701;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 3.28084;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 0.000621371;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 0.000539957;
            end
        case {"kilometer", "kilometers", "km"}
            switch desiredUnits
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 1000000;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 100000;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 1000;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 1;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 39370.1;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 3280.84;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 0.621371;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 0.539957;
            end
        case {"inch", "inches", "in"}
            switch desiredUnits
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 25.4;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 2.54;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 0.0254;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 2.54e-5;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 1;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 0.0833333;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 1.5783e-5;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 1.3715e-5;
            end
        case {"foot", "feet", "ft"}
            switch desiredUnits
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 304.8;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 30.48;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 0.3048;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 0.0003048;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 12;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 1;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 0.000189394;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 0.000164579;
            end
        case {"mile", "miles", "mi"}
            switch desiredUnits
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 160934000;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 160934;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 1609.34;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 1.60934;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 63360;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 5280;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 1;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 0.868976;
            end
        case {"nautical mile", "nautical miles", "nmi"}
            switch desiredUnits
                case {"millimeter", "millimeters", "mm"}
                    unitmult(ii) = 1852000;
                case {"centimeter", "centimeters", "cm"}
                    unitmult(ii) = 185200;
                case {"meter", "meters", "m"}
                    unitmult(ii) = 1852;
                case {"kilometer", "kilometers", "km"}
                    unitmult(ii) = 1.852;
                case {"inch", "inches", "in"}
                    unitmult(ii) = 72913.4;
                case {"foot", "feet", "ft"}
                    unitmult(ii) = 6076.12;
                case {"mile", "miles", "mi"}
                    unitmult(ii) = 1.15078;
                case {"nautical mile", "nautical miles", "nmi"}
                    unitmult(ii) = 1;
            end
    end
    
    %% Speed
    
    %% Mass
    
    %% Mass Flow
    
    %% Force
    
    %% Energy
    
end