function geopH = convertGeometricHeightToGeopotentialHeight(Reff, geomH)
% 
% Matt Werner (m.werner@vt.edu) - Dec 5, 2020
% 
% Calculate geopotential altitude based off a spherical Earth model of
% effective radius Reff at an altitude geomH. The spherical Earth has a
% gravitational field such that the gravitational acceleration at a
% geometric (physical) height geomH is
%                          /      Reff      \ 2
%            g(geomH) = g  | -------------- |   ,
%                        0 \  Reff + geomH  /
% where g0 is the gravitational acceleration at sea level at the effective
% radius of the Earth, Reff, which is the distance from the Earth's center
% to the point at sea level where the gravitational acceleration at sea
% level is precisely g = g0. The geopotential height is then calculated
% according to its definition.
%                          -1  / geomH
%               geopH =  g     |        g(geomH') dgeomH'.
%                         0    / 0
% 
%    Inputs:
% 
%              Reff - Effective radius of the Earth such that the
%                     gravitational acceleration experienced at sea-level
%                     is g = g0. The effective radius is simply the
%                     straight-line distance from the real Earth's center
%                     to the physical location at sea-level. Using Reff
%                     facilitates modeling Earth as a perfect sphere of
%                     radius Reff.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%             geomH - Geometric (physical) altitude above mean sea-level
%                     (MSL) to convert into geopotential altitude above
%                     mean sea-level.
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
% 
%    Outputs:
% 
%             geopH - Geopotential altitude above mean sea-level (MSL)
%                     calculated from the geometric altitude above MSL
%                     according to the Newtonian gravity model (perfect
%                     spherical Earth).
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
% 

%  Calculate the geopotential altitude
geopH = Reff*geomH / (Reff + geomH);
