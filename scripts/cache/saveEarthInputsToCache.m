%
% Matt Werner (m.werner@vt.edu) - Jan 9, 2021
% 
% SCRIPT - Create predetermined, temporary variables from the 'earth'
% structure to save into the cache as to avoid issues with the structure
% being dynamic and becoming (possibly) very large
% 

% Create temporary variables that are predetermined (hard-coded) for use
% within the cache to see if the requested Earth model for the current
% simulation has changed from the Earth model used for the previous
% simulation
earthModel = earth.model;
gravityModel = earth.gravity.model;
gravityDegree = earth.gravity.degree;
gravityOrder = earth.gravity.order;
magneticDegree = earth.magnetic.degree;
magneticOrder = earth.magnetic.order;
terrainAngleUnits = earth.terrain.angleUnits;

% Save the inputs that determine the Earth model to the cache
save(pathToCache, 'earthModel', 'gravityModel', 'gravityDegree', ...
    'gravityOrder', 'magneticDegree', 'magneticOrder', ...
    'terrainAngleUnits')

% Remove temporary variables
clear earthModel gravityModel gravityDegree gravityOrder ...
      magneticDegree magneticOrder terrainAngleUnits