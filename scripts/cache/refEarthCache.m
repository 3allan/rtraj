% 
% Matt Werner (m.werner@vt.edu) - Dec. 24, 2020
% 
% Cache system to improve the load time required to initialize a run of
% rtraj. The most significant portion of time, by far, is spent in loading
% the terrain.

% Set default flags indicating that no cache currently exists and
% all components of the model determining the ellipsoid parameters,
% gravitational and magnetic fields, and terriain shall be loaded as if for
% the first time
flags.exists.previous.EarthModelCache = false;
flags.exists.previous.rtrajWorkspace = true;
flags.load.ellipsoid = true;
flags.load.gravityField = true;
flags.load.magneticField = true;
flags.load.terrain = true;

% Check if the cache system should be used. If not, stop here; otherwise,
% continue
if (~flags.options.use.cache), return, end

% Attempt to load previous inputs for the earth model into the workspace.
% This process is justified because the load only occurs once per run of
% rtraj.m and the resulting table contains few entries, so the load is fast
try
    % Obtain the path that points to where the earth cache is
    pathToCache = getPathToCache('previous_earth_model.mat');
    % Attempt to load the cache (may not exist)
    previousEarth = load(pathToCache);
    flags.exists.previous.EarthModelCache = true;
catch error_CacheMiss
    % Check possible causes for error
    switch error_CacheMiss.identifier
        case 'MATLAB:load:couldNotReadFile'
            % Assume the workspace is empty since there is no cache
            flags.exists.previous.rtrajWorkspace = false;
            % Create a new cache file since one currently doesn't exist and
            % carry on to load the earth model
            saveEarthInputsToCache
        otherwise
            rethrow(error_CacheMiss)
    end
    return
end

% Check if the models exist in the rtraj.m workspace. To avoid issues with
% complexity and naming schemes, only check if the variable 'GM' exists
% and assume the others haven't been deleted from the previous run (if
% there was one). Therefore, assume that the other variables are there as
% well if just one of them is there since they all load together
try
    % Attempt to access the previous 'earth' structure's gravitational
    % parameter, GM
    earth.pars.GM; clear ans
    % Assume that nothing has been deleted from it since the last run, so
    % the previous workspace is functional as the cache itself
    flags.exists.previous.rtrajWorkspace = true;
catch
    % Regardless of whatever problem that occured, assume that the entire
    % workspace is gone. This way, the constituents of every model are
    % loaded. Don't throw an error
    flags.exists.previous.rtrajWorkspace = false;
end

% Only load the models that have changed from the previous run, if there
% was one
if (flags.exists.previous.rtrajWorkspace)
    % Determine which models require loading/processing by comparing currently
    % requested model parameters to previous model parameters
    % 
    % Check if the requested ellipsoid is different from the one used in the
    % previous run
    if (strcmp(earth.model, previousEarth.earthModel))
        flags.load.ellipsoid = false;
    end
    % Check if the requested gravitational field is different from the one
    % used in the previous run
    if (earth.gravity.degree == previousEarth.gravityDegree && ...
            earth.gravity.order == previousEarth.gravityOrder && ...
            strcmp(earth.gravity.model, previousEarth.gravityModel))
        flags.load.gravityField = false;
    end
    % Check if the requested magnetic field is different from the one used in
    % the previous run
    if (earth.magnetic.degree == previousEarth.magneticDegree && ...
            earth.magnetic.order == previousEarth.magneticOrder)
        flags.load.magneticField = false;
    end
    % Check if the requested angular units for the terrain are different from
    % the ones used in the previous run
    if (strcmp(earth.terrain.angleUnits, previousEarth.terrainAngleUnits))
        flags.load.terrain = false;
    end
end

% Update the cache
saveEarthInputsToCache