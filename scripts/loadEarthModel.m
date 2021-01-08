% 
% Matt Werner (m.werner@vt.edu) - Dec 28, 2020
% 
% Load the earth model while attempting to use the cache if possible.
% 

% SCRIPT - Creates flags to determine which models may be used again from
% the rtraj.m workspace (if already populated by a previous simulation) and
% which require reloading
refEarthCache

% Only load the models and coefficients used to define the earth that are
% either not currently stored in the rtraj.m workspace or are to be updated
% from the previous simulation
% 
% Load ellipsoid parameters
if (flags.load.ellipsoid)
    [GM, Req, Rpo, f, e, w] = defineEllipsoidParameters(earthModel);
end
% Load the geopotential harmonic coefficients
if (flags.load.gravityField)
    [nG, mG, ~, Cnm, Snm, ~, ~] = loadGravitationalCoefficients(gravityDegree, gravityOrder, gravityModel);
    [Cnm, Snm] = updateGravitationalCoefficients(JDLaunch, nG, mG, Cnm, Snm);
else
    disp("Hitting gravity model cache...")
    fprintf("Loaded gravity model\n\n")
end
% Load the magnetic field potential harmonic coefficients
if (flags.load.magneticField)
    [nM, mM, gnm, hnm, dgnmdt, dhnmdt] = loadMagneticCoefficients(magneticDegree, magneticOrder);
    [gnm, hnm] = updateMagneticCoefficients(JDLaunch, gnm, hnm, dgnmdt, dhnmdt);
else
    disp("Hitting magnetic model cache...")
    fprintf("Loaded magnetic model\n\n")
end
% Load the terrain
if (flags.load.terrain)
    [longitudes, geodeticLatitudes, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain(terrainAngleUnits);
else
    disp("Hitting terrain cache...")
    fprintf("Loaded terrain\n\n")
end

% Clear up workspace of temporary variables
clear previous_earth_model