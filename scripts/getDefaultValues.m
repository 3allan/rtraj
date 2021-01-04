% 
% Matt Werner (m.werner@vt.edu) - Jan 4, 2020
% 
% Provide default values and (relative) uncertainty for use within
% profile files using SI units.
% 
% SCRIPT - Adds ONE new variable to the current workspace called DEFAULT.
% 

% Define default values in SI units
% 
% Default relative (unitless) uncertainty (0.10 = 10%)
DEFAULT.RELATIVE_UNCERTAINTY = 0.10;
% 
% Default (atmospheric) ambient pressure (100000 Pa = 100 kPa)
DEFAULT.AMBIENT_PRESSURE = 100000; 
% 
% Default heat capacity ratio (ratio of specific heats) for the default gas
% (standard air)exhausting from the motor (1.40 = 7/5 = 1 + 2/(5 = #dof))
DEFAULT.EXHAUST_HEAT_CAPACITY_RATIO = 1.40;