%% Test how long it takes to get some variables
tic
for ii = 1:101
    defineEllipsoidParameters('WGS84');
end
toc

tic
for ii = 1:10000
    checkInput('WGS84');
end
toc

%% Test how long it takes to assign a variable
tic
testx = 1;
toc

tic
testx = rand(5000000,1);
toc

%% Test if it would be faster to implement jet atmosphere as exponential or polynomial
h = rand(5000, 1);

tic
exp(h);
toc

tic
p = 101325 - 12.0154*h + 0.000575298*h.^2 - 1.40015e-8*h.^3 + 1.76465e-13*h.^4;
toc

%% Test which method to obtain elevation above MSL is fastest so far (Dec 4, 2020)
clearvars
[longitude, geodeticLatitude, WGS84ToGeoid, GeoidToTerrain, WGS84ToTerrain] = loadTerrain("radians");
load("mats/EGM/coefficientmats/EGM2008HeightsCoeff_wrtGeoid.mat");
HN = EGM2008HeightsCoeff_wrtGeoid(:, 1);
HM = EGM2008HeightsCoeff_wrtGeoid(:, 2);
HC = EGM2008HeightsCoeff_wrtGeoid(:, 3);
HS = EGM2008HeightsCoeff_wrtGeoid(:, 4);
[~, Req, ~, f, ~, ~] = defineEllipsoidParameters('WGS84');

lonq = 4.3240;
latqgc = 15.21340; % Geocentric
latqgd = geoc2geod(latqgc, Req, 'WGS84');

tic
computeMSLtoTerrain(deg2rad(lonq), deg2rad(latqgc), HC, HS)
toc

tic 
computeSphericalHarmonicSeries(HN, HM, HC, HS, deg2rad(lonq), deg2rad(latqgc))
toc

tic
fastinterp2(longitude, geodeticLatitude, GeoidToTerrain, deg2rad(lonq), deg2rad(latqgd))
toc