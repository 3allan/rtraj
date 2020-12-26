% Cache system to improve the load time required to initialize a run of
% rtraj. The most significant portion of time, by far, is spent in loading
% the terrain.

% Create a flag to indicate that the loading process should load the
% terrain matrices, incurring lots of run time.
flag_loadTerrainMats = true;
if (exist('GeoidToTerrain', 'var'))
    % Switch the flag to be false
    falg_loadTerrainMats = false;
end