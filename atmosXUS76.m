function [d, T, p, c, mu] = atmosXUS76(geopH, gSL, Reff)
% 
% Matt Werner (m.werner@vt.edu) - Dec 6, 2020
% 
% Define the Extended 1976 U.S. Standard Atmosphere, but while allowing for
% a 'variable' effective radius of Earth, gravitational acceleration at sea
% level, and air gas constant.
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
%             Reff - Reference radius of the spherical model of Earth. This
%                    quantity does not affect the 1976 U.S. Standard
%                    Atmosphere in this function, but it does affect the
%                    Extended 1976 U.S. Standard Atmosphere. Therefore, it
%                    is unrequired if the geopotential altitude is beneath
%                    150 km altitude above mean sea level (MSL).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
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
    % Calculate density, temperature, and pressure according to US76.
    % density abides by the ideal gas law (p = d*R*T), pressure is a
    % decaying exponential, and temperature is a (piecewise continuous)
    % linear function of geopotential altitude. The temperature that is
    % output from this atmosphere model is actually molecular temperature,
    % but for altitudes still within the majority of the atmosphere, the
    % molecular and environmental temperatures are the same.
    
    % Define gas constant of dry air as R = 287.053 (atmos.m [03], table 4)
    
    % Troposphere
    if (geopH < 11000)
    % (geopH_, d_, T_, p_, dTdgeopH)=(0, 1.225, 288.15, 101325, -0.0065)
        [d, T, p] = computeSAP(geopH, 0, 1.225, 288.15, 101325, -0.0065, gSL, 287.053);
    % Stratosphere
    elseif (geopH < 20000)
        % (geopH_, d_, T_, p_, dTdgeopH)=(11000, 0.363717160273004, 216.65, 22619.582, 0)
        [d, T, p] = computeSAP(geopH, 11000, 0.363717160273004, 216.65, 22619.582, 0, gSL, 287.053);
    elseif (geopH < 32000)
        % (geopH_, d_, T_, p_, dTdgeopH)=(20000, 0.0879639515286987, 216.65, 5470.481, 0.001)
        [d, T, p] = computeSAP(geopH, 20000, 0.0879639515286987, 216.65, 5470.481, 0.001, gSL, 287.053);
    elseif (geopH < 47000)
        % (geopH_, d_, T_, p_, dTdgeopH)=(32000, 0.0132142969635127, 228.65, 867.316, 0.0028)
        [d, T, p] = computeSAP(geopH, 32000, 0.0132142969635127, 228.65, 867.316, 0.0028, gSL, 287.053);
    % Mesosphere
    elseif (geopH < 51000)
        % (geopH_, d_, T_, p_, dTdgeopH)=(47000, 0.001428610145188, 270.65, 110.990, 0)
        [d, T, p] = computeSAP(geopH, 47000, 0.001428610145188, 270.65, 110.990, 0, gSL, 287.053);
    elseif (geopH < 71000)
        % (geopH_, d_, T_, p_, dTdgeopH)=(51000, 0.000862211723448629, 270.65, 66.986, -0.0028)
        [d, T, p] = computeSAP(geopH, 51000, 0.000862211723448629, 270.65, 66.986, -0.0028, gSL, 287.053);
    else % so geopH < 84852 m,
        % (geopH_, d_, T_, p_, dTdgeopH)=(71000, 6.37822106462866e-05, 214.65, 3.930, -0.0020)
        [d, T, p] = computeSAP(geopH, 71000, 6.37822106462866e-05, 214.65, 3.930, -0.0020, gSL, 287.053);
    end
    % Stop before the thermosphere
else
    % Calculate the Extended 1976 U.S. Standard Atmosphere 
    geomHkm = 0.001*geomH;
    Reffkm = 0.001*Reff;
    geomHkmPowers = [geomHkm^4, geomHkm^3, geomHkm^2, geomHkm, 1];
    % Thermosphere
    if (geomHkm < 91)
        P = [0.000000, 	2.159582E-06, 	-4.836957E-04, 	-0.1425192, 	13.47530]';
        R = [0.000000, 	-3.322622E-06, 	9.111460E-04, 	-0.2609971, 	5.944694]';
        T = 186.8673;
    elseif (geomHkm < 100)
        P = [0.000000, 	3.304895E-05, 	-0.009062730, 	0.6516698, 	-11.03037]';
        R = [0.000000, 	2.873405E-05, 	-0.008492037, 	0.6541179, 	-23.62010]';
        T = 263.1905 - 76.3232 * sqrt(1 - ((geomHkm - 91) / -19.9429)^2);
    elseif (geomHkm < 110)
        P = [0.000000, 	6.693926E-05, 	-0.01945388, 	1.719080, 	-47.75030]';
        R = [-1.240774E-05, 	0.005162063, 	-0.8048342, 	55.55996, 	-1443.338]';
        T = 263.1905 - 76.3232 * sqrt(1 - ((geomHkm - 91) / -19.9429)^2);
    elseif (geomHkm < 120)
        P = [0.000000, 	-6.539316E-05, 	0.02485568, 	-3.223620, 	135.9355]';
        R = [0.00000, 	-8.854164E-05, 	0.03373254, 	-4.390837, 	176.5294]';
        T = 240 + 12 * (geomHkm - 110);
    elseif (geomHkm < 150)
        P = [2.283506E-07, 	-1.343221E-04, 	0.02999016, 	-3.055446, 	113.5764]';
        R = [3.661771E-07, 	-2.154344E-04, 	0.04809214, 	-4.884744, 	172.3597]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
    elseif (geomHkm < 200)
        P = [1.209434E-08, 	-9.692458E-06, 	0.003002041, 	-0.4523015, 	19.19151]';
        R = [1.906032E-08,	-1.527799E-05, 	0.004724294, 	-0.6992340, 	20.50921]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
    elseif (geomHkm < 300)
        P = [8.113942E-10, 	-9.822568E-07, 	4.687616E-04, 	-0.1231710, 	3.067409]';
        R = [1.199282E-09, 	-1.451051E-06, 	6.910474E-04, 	-0.1736220, 	-5.321644]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
    elseif (geomHkm < 500)
        P = [9.814674E-11, 	-1.654439E-07, 	1.148115E-04, 	-0.05431334, 	-2.011365]';
        R = [1.140564E-10, 	-2.130756E-07, 	1.570762E-04, 	-0.07029296, 	-12.89844]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
        % ~Thermosphere/Exosphere
    elseif (geomHkm <= 750)
        P = [-7.835161E-11, 	1.964589E-07, 	-1.657213E-04, 	0.04305869,      -14.77132]';
        R = [8.105631E-12, 	-2.358417E-09, 	-2.635110E-06, 	-0.01562608, 	-20.02246]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
        % Exosphere
    elseif (geomHkm < 1000)
        P = [2.813255E-11, 	-1.120689E-07, 	1.695568E-04, 	-0.1188941,      14.56718]';
        R = [-3.701195E-12, 	-8.608611E-09, 	5.118829E-05, 	-0.06600998, 	-6.137674]';
        xi_expfactor = (geomHkm - 120) * (Reffkm + 120) / (Reffkm + geomHkm);
        T = 1000 - 640*exp(-0.01875*xi_expfactor);
    else
        % Beyond model
        error("No atmospheric data for altitudes exceeding 1000 km (Extended 1976 USSA).")
    end
    % Calculate exponential power series of density and pressure
    d = exp(geomHkmPowers*R);
    p = exp(geomHkmPowers*P);
end


% Compute the speed of sound (sqrt(1.4*287.063) = 20.0471494232971)
c = computeSOS(1.4, 287.053, T);

% Compute the dynamic viscosity according to Sutherland's law
mu = computeDynamicViscosity(T);