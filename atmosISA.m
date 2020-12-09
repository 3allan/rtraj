function [d, T, p, c, mu] = atmosISA(geopH, gSL)
% 
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Define the International Standard Atmosphere (ISA) up to a geometric
% (physical) altitude of 86,000 meters above mean sea level (MSL). Given
% the (base) pressure, temperature, and density at the beginning of each
% atmospheric layer from -610 meters below MSL to 86,000 meters above MSL,
% the specific gas constant of air can be found using the ideal gas law and
% varied between atmospheric layers. The gas constant stays relatively
% constant (~287 J/kgK) from sea level until reaching the mesosphere, where
% it drops (~205 J/kgK). Above this point, the density is absent from the
% presented base values, so it assumed to remain constant (~205 J/kgK) from
% this point onwards until meeting the 86,000 meter ceiling. The gas
% constant is given a linear profile in geopotential altitude from sea
% level until 86,000 according to the rule of mixtures, piecing adjacent
% gas constants together (piecewise-) linearly.
% 
% 
%    Inputs:
% 
%             geopH - Geopotential height above mean sea level (MSL). Given
%                     the geopotential height, the temperature profile (in
%                     a particular layer) of the standard atmosphere can be
%                     calculated as a simple linear variation, depending
%                     only on the atmospheric temperature at the base of
%                     the layer, the lapse rate, and the difference between
%                     geopotential altitudes of the current position and
%                     the beginning of the layer.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%               gSL - Gravitational acceleration (positive) at sea level
%                     directly beneath the current position. For standard
%                     atmospheres that take Earth to be a perfect sphere,
%                     this value is constant and, for the 1976 USSA model,
%                     is taken to be exactly 9.80665. Otherwise, this
%                     quantity may vary with geodetic latitude (as would
%                     the geopotential altitude) if attempting to use the
%                     standard atmosphere numbers with an ellipsoidal
%                     Earth.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s2 (meters per squared second)
% 
%    Outputs:
% 
%                 d - Atmospheric density.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
%                 T - Environmental (kinetic) temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%                 p - Atmospheric pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%         varargout - Additional outputs specified.
%                     Size: ?
%                     Units: ?
% 

% Calculate atmospheric properties at the given geopotential altitude
if (geopH < 84852)
    % Calculate density, temperature, and pressure according to ISA.
    % density abides by the ideal gas law (p = d*R*T), pressure is a
    % decaying exponential, and temperature is a (piecewise continuous)
    % linear function of geopotential altitude. The temperature that is
    % output from this atmosphere model is actually molecular temperature,
    % but for altitudes still within the majority of the atmosphere, the
    % molecular and environmental temperatures are basically the same.
    
    % Troposphere
    if (geopH < 11000)
        Rair = mixf((geopH + 610)/11610, [287.065, 285.747]);
        % (geopH_, d_, T_, p_, dTdgeopH)=(-610, 1.2985, 292.15, 108900, -0.0065)
        [d, T, p] = computeSAP(geopH, -610, 1.2985, 292.15, 108900, -0.0065, gSL, Rair);
    % Tropopause
    elseif (geopH < 20000)
        Rair = mixf((geopH - 11000)/9000, [285.747, 285.848]);
        % (geopH_, d_, T_, p_, dTdgeopH)=(11000, 0.3639, 217.65, 22632, 0)
        [d, T, p] = computeSAP(geopH, 11000, 0.3639, 217.65, 22632, 0, gSL, Rair);
    elseif (geopH < 32000)
        Rair = mixf((geopH - 20000)/12000, [285.848, 287.597]);
        % (geopH_, d_, T_, p_, dTdgeopH)=(20000, 0.0880, 217.65, 5474.9, 0.001)
        [d, T, p] = computeSAP(geopH, 20000, 0.0880, 217.65, 5474.9, 0.001, gSL, Rair);
    elseif (geopH < 47000)
        Rair = mixf((geopH - 32000)/15000, [287.597, 204.9]);
        % (geopH_, d_, T_, p_, dTdgeopH)=(32000, 0.0132, 228.65, 868.02, 0.0028)
        [d, T, p] = computeSAP(geopH, 32000, 0.0132, 228.65, 868.02, 0.0028, gSL, Rair);
    % Mesosphere
    elseif (geopH < 51000)
        Rair = 204.9;
        % (geopH_, d_, T_, p_, dTdgeopH)=(47000, 0.0020, 270.65, 110.91, 0)
        [d, T, p] = computeSAP(geopH, 47000, 0.0020, 270.65, 110.91, 0, gSL, Rair);
    elseif (geopH < 71000)
        Rair = 204.9;
        % (geopH_, d_, T_, p_, dTdgeopH)=(51000, 0.0012, 270.65, 66.939, -0.0028)
        [d, T, p] = computeSAP(geopH, 51000, 0.0012, 270.65, 66.939, -0.0028, gSL, Rair);
    else % so geopH < 84.852,
        Rair = 204.9;
        % (geopH_, d_, T_, p_, dTdgeopH)=(71000, 0.000089955, 214.65, 3.9564, -0.0020)
        [d, T, p] = computeSAP(geopH, 71000, 0.000089955, 214.65, 3.9564, -0.0020, gSL, Rair);
    end
    % Stop before the thermosphere
else
    % Beyond model
    error("No atmospheric data for altitudes exceeding 1000 km (ISA).")
end

% Compute the speed of sound (sqrt(1.4*287.063) = 20.0471494232971)
c = computeSOS(1.4, Rair, T);

% Compute the dynamic viscosity according to Sutherland's law
mu = computeDynamicViscosity(T);