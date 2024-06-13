% Generate orthometric heights H as specified by EGM2008


geocentricLatitude = linspace(-pi/2, pi/2, 4321)';
longitude = linspace(0, 2*pi*8639/8640, 8640);

load("coefficientmats/EGM2008HeightsCoeff_wrtGeoid.mat", "EGM2008HeightsCoeff_wrtGeoid");
% n = EGM2008HeightsCoeff_wrtGeoid(:, 1);
% m = EGM2008HeightsCoeff_wrtGeoid(:, 2);
HC = EGM2008HeightsCoeff_wrtGeoid(:, 3);
HS = EGM2008HeightsCoeff_wrtGeoid(:, 4);

cd("../../")
tic
H = computeMSLtoTerrain(HC, HS, longitude, geocentricLatitude);
toc

% SAVE THE FILE
1
1
1
1