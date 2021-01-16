%% Test obtaining the nose cone properties (again, new function)
clear, clc
nosecone.OD = 2*0.0254*3;
nosecone.length = 0.0254*30;
nosecone.thickness = 0.0254*0.25;
nosecone.k = 0;
nosecone.density = 2700;
nosecone.series = "haack";
nosecone.thicknessType = "vertical";

[M, CoM, IMoIatCoM, profile] = ...
    computeNoseconeProperties(nosecone);

x = profile(:, 1);
outersurf = profile(:, 2);
innersurf = profile(:, 3);

plot(x, outersurf, 'b', x, -outersurf, 'b')
hold on
tmpinnersurf = innersurf;
tmpinnersurf(innersurf == 0) = NaN;
plot(x, tmpinnersurf, 'r', x, -tmpinnersurf, 'r')
hold off
axis equal

X = [flip(x); x];
Y = [flip(outersurf); -outersurf];
YY = [flip(innersurf); -innersurf];
fill(X, Y, 'b')
hold on
fill(X, YY, 'w')
hold off
axis equal


%% Test cosine spacing
clear, clc
L=1;
L=5;
x = linspace(0, L, 10000);
t = linspace(0, pi, 10000);
xx = 0.5*L*(1 - cos(t));
plot(x, xx)

%% Test BOR properties and ensure they match well with results from below

R = 0.0762;
L = 0.762;
l = L;
k = 0;

x = linspace(0, L, 500000)';
f = R*(x/L) .* (2 - k*(x/L))/(2 - k);
g = R*((x-l)/L) .* (2 - k*((x-l)/L))/(2 - k);
g(x <= l) = 0;

plot(x, f, x, g)

[M, CoM, IMoIatCoM] = computeBORProperties(x, f, g, 2700)

1;

%% Test time it takes to calculate structural properties & see if they're correct
% R = 1;
% L = 5*R;
% l = R/40;
% k = 0.845;

R = 0.0762;
L = 0.762;
l = L;
k = 0;

tic
V = computeParabolicConeVolume(R, L, l, k)
M = V*2700
CoM = computeParabolicConeCenterOfMass(L, l, k)
I = computeParabolicConeInertiaMatrix(R, L, l, k, M, CoM)
toc

% l = 0 matches solid cone

%% Test if algorithm for determining which profiles receive Monte Carlo treatment is correct
A = randi(2, 2, 4)-1;
B = randi(2, 2, 5)-1;
C = B - A;

idxNoDiffs = find(C == 0);
C(idxNoDiffs) = C(idxNoDiffs) + A(idxNoDiffs);

% Change all -1 values to 0 to indicate that the values are off for
% this simulation
C(C == -1) = 0;

flag = any(C, [1,2]);


%% Test string equivalence of some string arrays

% Should pass
testStr = ["aaaa", "aaaa", "aaaa", "aaaa"];
checkStringEquivalence(testStr)

% Should fail at position (1,2)
% testStr = ["aaaa", "abaa", "aaaa", "aaaa"];
% checkStringEquivalence(testStr)

% Should fail at position (3,2)
testStr = ["aaaa", "aaaa", "aaOa", "aaaa" 
           "aaaa", "aaaa", "aaaa", "aaaa"
           "aaaa", "aOaa", "aaaa", "aaaa"
           "aaaa", "aaaa", "aOaa", "aaaa"
           "aaaa", "aaaa", "aaaa", "aaaa"];
checkStringEquivalence(testStr)


%% Test if it's faster to check a condition than assign large variable
tic
if (1 > 2)
end
toc

tic
testx = rand(100, 1);
toc


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