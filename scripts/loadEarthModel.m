% 
% Matt Werner (m.werner@vt.edu) - Dec 28, 2020
% 
% Load the earth model while attempting to use the cache if possible.

% SCRIPT - Creates flags to determine which models may be used again from
% the rtraj.m workspace (if already populated by a previous simulation) and
% which require reloading
cache

% Only load the models and coefficients usd to define the earth that are
% either not currently stored in the rtraj.m workspace or are to be updated
% from the previous simulation
% 
% Load or cache ellipsoid parameters
if (flag_loadEllipsoid)
    [GM, Req, Rpo, f, e, w] = defineEllipsoidParameters(earthModel);
end
% Load or cache the geopotential harmonic coefficients
if (flag_loadGravityField)
    [nG, mG, ~, Cnm, Snm, ~, ~] = loadGravitationalCoefficients(gravityDegree, gravityOrder, gravityModel);
    [Cnm, Snm] = updateGravitationalCoefficients(JDLaunch, nG, mG, Cnm, Snm);
else
    disp("Caching gravity model...")
    fprintf("Loaded gravity model\n\n")
end
% Load or cache the magnetic field potential harmonic coefficients
if (flag_loadMagneticField)
    [nM, mM, gnm, hnm, dgnmdt, dhnmdt] = loadMagneticCoefficients(magneticDegree, magneticOrder);
    [gnm, hnm] = updateMagneticCoefficients(JDLaunch, gnm, hnm, dgnmdt, dhnmdt);
else
    disp("Caching magnetic model...")
    fprintf("Loaded magnetic model\n\n")
end
% Load or cache the terrain
if (flag_loadTerrain)
    [longitudes, geodeticLatitudes, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain(terrainAngleUnits);
else
    disp("Caching terrain...")
    fprintf("Loaded terrain\n\n")
end