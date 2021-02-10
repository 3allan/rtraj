% 
% Feb. 10, 2021
% 
% Define the parameters for determining an Earth model as a Matlab
% structure that holds information about the Earth's reference ellipsoid,
% gravitational and magnetic field models, and conditions regarding the
% launch site (time and location, etc.)
% 

% Obtain actual time right now
earth.time.now = datetime('now', 'timezone', 'America/New_York');
% Obtain actual launch times as Julian date (UT1) and other timezones
% (UTC & local)
[earth.time.launch.JD, ...
 earth.time.launch.UTC, ...
 earth.time.launch.LOC] ...
 = defineLaunchTime('15-July-2032 13:53:45.234', 'd-MMMM-yyyy HH:mm:ss.SSS', '-05:00');

% Provide earth model parameters
earth.model = "WGS84"; % Ellipsoid
earth.gravity.model = "tide-free"; % tide-free or zero-tide (EGM08)
earth.gravity.degree = 15; % Degree of gravity model (EGM08)
earth.gravity.order = 15; % Order of gravity model (EGM08)
earth.magnetic.degree = 12; % Degree of magnetic model (WMM20)
earth.magnetic.order = 12; % Order of magnetic model (WMM20)
earth.atmos.model = "NRLMSISE00"; % Atmosphere model
earth.terrain.angleUnits = "degrees"; % Units for lat/lon on terrain map (plotting only)
% Provide launch site location
earth.launchsite.longitude = convUnits(-119.0560102, "degrees", "radians"); % Longitude of launch site [rad]
earth.launchsite.latitude  = convUnits(+ 40.9107330, "degrees", "radians"); % Geodetic latitude of launch site [rad]
% Provide launch site temperature and density
earth.launchsite.temperature = 70; % Local surface temperature [F]
% Define launch tower length
earth.launchsite.towerLength = convUnits(30, "feet", "meters"); % Length of launch tower [m]