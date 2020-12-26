function dynamicViscosity = computeDynamicViscosity(T)
% 
% Matt Werner (m.werner@vt.edu) - Dec 5, 2020
% 
% Calculate dynamic viscosity (standard symbol - mu) of dry, ideal air at a
% temperature T according to Sutherland's law,
%                             3/2
%                   /   T    \      Tref + S
%       mu = muref  | ------ |     ----------,
%                   \  Tref  /       T + S
% where muref = 0.00001716 Pa*s (Pascal seconds) and Tref = 273.15 K
% (Kelvin) are the reference viscosity and temperature values and S = 110.4
% K (Kelvin) is Sutherland's constant.
% 
%    Inputs:
% 
%                 T - Air temperature.
%                     Size: n-by-1 (vector)
%                     Units: K (Kelvin)
% 
%    Outputs:
% 
%  dynamicViscosity - Dynamic viscosity of air at temperature T.
%                     Size: n-by-1 (vector)
%                     Units: Pa*s (Pascal seconds)
% 

% Calculate dynamic viscosity as based off of Sutherland's law
dynamicViscosity = 0.00000145793265 * T.^1.5 ./ (T + 110.4);
