function [d, T, p] = computeSAP(geopH, geopH_, d_, T_, p_, dTdgeopH, gSL, R)
%
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Calculate the (S)tandard (A)tmosphere (P)roperties of any particular
% standard atmosphere model defined by the base properties d_, T_, and p_
% of a layer of the atmosphere at the geopotential height geopH_. The
% gravitational acceleration at sea level is g0 and the constituents of air
% are encapsulated by the air's gas constant R.
% 
%    Inputs:
% 
%             geopH - Geopotential altitude above the effective radius of
%                     Earth, be it constant (1976 USSA, for example) or
%                     varying with latitude.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%            geopH_ - Geopotential altitude of the beginning of the current
%                     atmospheric layer which is determined by the physical
%                     altitude above sea level.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                d_ - Initial density at the lowest level of the current
%                     atmospheric layer at geopotential height. The initial
%                     value at sea level is usually ~1.225 kg/m3 (1976
%                     USSA).
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%                T_ - Initial temperature at the lowest level of the
%                     current atmospheric layer at geopotential height. The
%                     initial value at sea level is usually ~288.15 K (1976
%                     USSA).
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%                p_ - Initial pressure at the lowest level of the current
%                     atmospheric layer at geopotential height. The initial
%                     value at sea level is usually ~1.225 kg/m3 (1976
%                     USSA).
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%          dTdgeopH - Lapse rate of temperature with respect to
%                     geopotential altitude. If zero, then the atmospheric
%                     layer is said to be isothermal.
%                     Size: 1-by-1 (scalar)
%                     Units: K/m (Kelvin per meter)
% 
%               gSL - Gravitational acceleration (positive) at sea level.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per squared second)
% 
%                 R - Gas constant of air.
%                     Size: 1-by-1 (scalar)
%                     Units: J/kg K (Joules per (kilogram x Kelvin))
% 

% Precalculate the difference between the current geopotential altitude and
% the geopotential altitude at which the current atmospheric layer started
deltageopH = geopH - geopH_;

% Check if this difference is negative (zero or negative)
if (deltageopH == 0)
    d = d_;
    T = T_;
    p = p_;
    return
elseif (deltageopH < 0)
    error("Mismatched atmospheric layers (deltageopH = %1.2f m)", deltageopH);
end

% Check if the layer is isothermal or not
if (dTdgeopH == 0)
    T = T_;
    pOverp_ = exp(-gSL * deltageopH / (R*T));
    dOverd_ = pOverp_;
else
    T = T_ + dTdgeopH*deltageopH;
    TOverT_ = T/T_;
    pOverp_ = TOverT_^-(gSL / (dTdgeopH * R));
    dOverd_ = pOverp_/TOverT_;
end

% Calculate pressure and density
p = pOverp_*p_;
d = dOverd_*d_;
    
