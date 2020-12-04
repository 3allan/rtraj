function [ft2in, fps2mph, in2m, lbf2N, lb2kg, psi2Pa] = getUnitConversions
% 
% Matt Werner (m.werner@vt.edu) - Dec 1, 2020
% 
% Obtain common unit conversions between [US] to [US] and [US] to [SI] unit
% systems. The inverse relations are obtained upon simple inversion.
ft2in = 12; % [ft] --> [in]
fps2mph = 0.681818; % [ft/s] --> [mi/hr]
in2m = 0.0254; % [in] --> [m]
lbf2N = 4.44822; % [lbf] --> [N]
lb2kg = 0.453592; % [lb] --> [kg]
psi2Pa = 6894.76; % [psi] --> [Pa]
end