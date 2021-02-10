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
    [earth.pars.GM, ...
      earth.pars.Req, ...
      earth.pars.Rpo, ...
      earth.pars.f, ...
      earth.pars.e, ...
      earth.pars.w] = defineEllipsoidParameters(earth.model);
end
% Load the geopotential harmonic coefficients
if (flags.load.gravityField)
    tmp_JD = earth.time.launch.JD;
    [n, m, ~, Cnm, Snm, ~, ~] = ...
        loadGravitationalCoefficients(earth.gravity.degree, ...
                            earth.gravity.order, earth.gravity.model);
    [Cnm, Snm] = updateGravitationalCoefficients(tmp_JD, n, m, Cnm, Snm);
    clear tmp_JD
    % Add degrees/orders and the coefficients to the earth structure and
    % immediately remove them from the workspace
    earth.gravity.coefficients.degrees = n; clear n
    earth.gravity.coefficients.orders = m; clear m
    earth.gravity.coefficients.cosine = Cnm; clear Cnm
    earth.gravity.coefficients.sine = Snm; clear Snm
else
    disp("Hitting gravity model cache...")
    fprintf("Loaded gravity model\n\n")
end
% Load the magnetic field potential harmonic coefficients
if (flags.load.magneticField)
    tmp_JD = earth.time.launch.JD;
    [n, m, gnm, hnm, dgnmdt, dhnmdt] = ...
        loadMagneticCoefficients(earth.magnetic.degree, earth.magnetic.order);
    [gnm, hnm] = updateMagneticCoefficients(tmp_JD, gnm, hnm, dgnmdt, dhnmdt);
    clear tmp_JD dgnmdt dhnmdt
    % Add degrees/orders and the coefficients to the magnetic structure and
    % immediately remove them from the workspace
    earth.magnetic.coefficients.degrees = n; clear n
    earth.magnetic.coefficients.orders = m; clear m
    earth.magnetic.coefficients.cosine = gnm; clear gnm
    earth.magnetic.coefficients.sine = hnm; clear hnm
else
    disp("Hitting magnetic model cache...")
    fprintf("Loaded magnetic model\n\n")
end
% Load the terrain
if (flags.load.terrain)
    [earth.terrain.longitudes, earth.terrain.geodeticLatitudes, ...
        earth.terrain.WGS84ToGeoid, earth.terrain.GeoidToTerrain, ...
        earth.terrain.WGS84ToTerrain] ...
        = loadTerrain("radians");
else
    disp("Hitting terrain cache...")
    fprintf("Loaded terrain\n\n")
end