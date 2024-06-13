f = 1/298.257223563; % WGS84 Flattening []

load("mats/EGM/GeoidTopographyUndulation.mat", "Geoid_Topography_height", "longitude", "geocentricLatitude", "-mat");
Geoid2Topography = Geoid_Topography_height; % 
geocentricLongitude02piG2T = longitude; % longitude
geocentricColatitudeG2T = geocentricLatitude; %CO-latitude
clear Geoid_Topography_height longitude


load("mats/EGM/WGS84GeoidUndulation.mat", "WGS84_Geoid_h", "longitude", "geodeticLatitude", "-mat");
WGS842Geoid = WGS84_Geoid_h;
geodeticLongitudeW2G = longitude;
geodeticLatitudeW2G = geodeticLatitude;
clear WGS84_Geoid_h longitude

% interp2(geodeticLongitude, geodeticLatitude, WGS842Geoid, 0, 0)

% Set geocentric colatitude (0, pi) to (-pi/2, pi/2) by subtracting it from pi/2
geocentricLatitudeG2T = pi/2 - geocentricColatitudeG2T;
% Convert geocentric latitude to geodetic latitude at zero height (PVR Rn is
% unnecessary with zero height from the ellipsoid).
geodeticLatitudeG2T = atan2(sin(geocentricLatitudeG2T), cos(geocentricLatitudeG2T)*(1 - f)^2);
clear geocentricColatitudeG2T geocentricLatitudeG2T
% Change longitude components so that longitude lies within (-pi, pi) with
% 0 situated directly in the center
geocentricLongitude0pix2G2T = mod(geocentricLongitude02piG2T, pi);
% Find the index in which the set resets back to 0 - starting here to the
% end is the west side of the world.
numelLongG2T = numel(geocentricLongitude0pix2G2T);
% No checks needed - guaranteed to reset since long. was given on (0, 2pi)
for ii = 1:numelLongG2T
    if (geocentricLongitude0pix2G2T(ii + 1) < geocentricLongitude0pix2G2T(ii))
        idx_resetTo0 = ii + 1;
        break
    end
end
% Negate and flip the longitudes over to get western longitudes (0, pi)
geocentricLongitudeWestpi0G2T = -flip(geocentricLongitude0pix2G2T(idx_resetTo0+1:end));
LongitudeG2T = [geocentricLongitudeWestpi0G2T, geocentricLongitude0pix2G2T(1:idx_resetTo0)];
LongitudeG2T(end) = pi;
% Repeat flips with the full Z matrix
WestG2T = Geoid2Topography(:, idx_resetTo0+1:end); % Flip left/right
G2T = [WestG2T, Geoid2Topography(:, 1:idx_resetTo0)];













% Repeat process for W2G
clear idx_resetTo0
% Set geodetic latitude (-90, 90) to (-pi/2, pi/2)
geodeticLatitudeW2G = deg2rad(geodeticLatitudeW2G);
% Change longitude components so that longitude lies within (-pi, pi) with
% 0 situated directly in the center
geodeticLongitude0pix2W2G = mod(deg2rad(geodeticLongitudeW2G), pi);
% Find the index in which the set resets back to 0 - starting here to the
% end is the west side of the world.
numelLongW2G = numel(geodeticLongitude0pix2W2G);
% No checks needed - guaranteed to reset since long. was given on (0, 2pi)
for ii = 1:numelLongW2G
    if (geodeticLongitude0pix2W2G(ii + 1) < geodeticLongitude0pix2W2G(ii))
        idx_resetTo0 = ii + 1;
        break
    end
end
% Negate and flip the longitudes over to get western longitudes (0, pi)
geodeticLongitudeWestpi0W2G = -flip(geodeticLongitude0pix2W2G(idx_resetTo0+1:end));
LongitudeW2G = [geodeticLongitudeWestpi0W2G, geodeticLongitude0pix2W2G(1:idx_resetTo0)];
LongitudeW2G(end) = pi;
% Repeat flips with the full Z matrix
WestW2G = WGS842Geoid(:, idx_resetTo0+1:end); % Flip left/right
W2G = [WestW2G, WGS842Geoid(:, 1:idx_resetTo0)];













% Now interpolate both to have the same long and lats and add those two
% matrices together to obtain WGS84 to topography
interpLat = geodeticLatitudeG2T;
interpLon = LongitudeG2T;
[X, Y] = meshgrid(LongitudeW2G, geodeticLatitudeW2G);
[Xq, Yq] = meshgrid(interpLon, interpLat);
W2GwithInterpolatedG2Tcoords = interp2(X, Y, W2G, Xq, Yq, 'makima');


W2T = W2GwithInterpolatedG2Tcoords + G2T;
image(LongitudeG2T, geodeticLatitudeG2T, W2T);


%% Save new mat files
% Longitude_WGS84ToGeoid = LongitudeW2G;
% GeodeticLatitude_WGS84ToGeoid = geodeticLatitudeW2G;
% WGS84ToGeoid = W2G;
return
save("wgs84-to-egm08geoid_undulation=meters_LatLon=radians.mat", "LongitudeW2G", "geodeticLatitudeW2G", "W2G", '-v7')
LongitudeW2G = rad2deg(LongitudeW2G);
geodeticLatitudeW2G = rad2deg(geodeticLatitudeW2G);
save("wgs84-to-egm08geoid_undulation=meters_LatLon=degrees.mat", "LongitudeW2G", "geodeticLatitudeW2G", "W2G", '-v7')




save("egm08geoid-to-ground_orthometricHeight=meters_LatLon=radians.mat", "LongitudeG2T", "geodeticLatitudeG2T", "G2T", '-v7')
LongitudeG2T = rad2deg(LongitudeG2T);
geodeticLatitudeG2T = rad2deg(geodeticLatitudeG2T);
save("egm08geoid-to-ground_orthometricHeight=meters_LatLon=degrees.mat", "LongitudeG2T", "geodeticLatitudeG2T", "G2T", '-v7')





LongitudeW2GwithInterpolatedG2Tcoords = interpLon;
geodeticLatitudeW2GwithInterpolatedG2Tcoords = interpLat;
save("wgs84-to-egm08geoid_undulation=meters_LatLon=(egm08geoid-to-ground_LatLon=radians).mat", "LongitudeW2GwithInterpolatedG2Tcoords", "geodeticLatitudeW2GwithInterpolatedG2Tcoords", "W2GwithInterpolatedG2Tcoords", '-v7')
LongitudeW2GwithInterpolatedG2Tcoords = rad2deg(LongitudeW2GwithInterpolatedG2Tcoords);
geodeticLatitudeW2GwithInterpolatedG2Tcoords = rad2deg(geodeticLatitudeW2GwithInterpolatedG2Tcoords);
save("wgs84-to-egm08geoid_undulation=meters_LatLon=(egm08geoid-to-ground_LatLon=degrees).mat", "LongitudeW2GwithInterpolatedG2Tcoords", "geodeticLatitudeW2GwithInterpolatedG2Tcoords", "W2GwithInterpolatedG2Tcoords", '-v7')




LongitudeW2T = geodeticLongitudeW2G;
geodeticLatitudeW2T = geodeticLatitudeW2G;
save("wgs84-to-ground_ellipsoidalHeight=meters_LatLon=radians.mat", "LongitudeW2T", "geodeticLatitudeW2T", "W2T", '-v7')
LongitudeW2T = rad2deg(LongitudeW2T);
geodeticLatitudeW2T = rad2deg(geodeticLatitudeW2T);
save("wgs84-to-ground_ellipsoidalHeight=meters_LatLon=degrees.mat", "LongitudeW2T", "geodeticLatitudeW2T", "W2T", '-v7')


