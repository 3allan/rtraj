function c = computeSOS(Y, R, T)
% 
% Matt Werner (m.werner@vt.edu) - Dec 8, 2020
% 
% Calculate the (S)peed (o)f (S)ound according to the ideal gas law,
%                               ______
%                         c = \/Y R T ,
% where Y is the heat capacity ratio (ratio of specific heats), R is the
% gas constant specific to the medium through which sound is travelling,
% and T is the molecular temperature of the gas.
% 
%    Inputs:
% 
%                 Y - Heat capacity ratio.
%                     Size: n-by-1 (scalar)
%                     Units: - (unitless)
% 
%                 R - Gas constant of the particular medium through which
%                     sound is travelling.
%                     Size: n-by-1 (scalar)
%                     Units: J/kg K (Joules per (kilogram x Kelvin))
% 
%                 T - Molecular temperature of the gas.
%                     Size: n-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%    Outputs:
% 
%                 c - Local speed of sound according to the ideal gas law,
%                     dependent solely upon the heat capacity ratio Y, gas
%                     constant R, and molecular temperature T.
%                     Size: n-by-1 (scalar)
%                     Units: m/s (meters per second)
% 

% No checks necessary

% Calculate the speed of sound
c = sqrt(Y.*R.*T);