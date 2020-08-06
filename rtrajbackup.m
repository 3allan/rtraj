% OLVT - Trajectory and Analysis
% Simulate trajectory of rocket's uppermost stage
% Created - Oct 16, 2019
% v5.1.5 Apr 18, 2020
% Matt Werner (mwerner4@vt.edu)
clear variables; clc; close all

% Notes
% 1. Use a mono font to correctly display formatting of comments/tables
% 2. Disable wrap-lines in Command Window preferences to view tables easily
% 3. No support to integrate automatically on an interval [t0, t1] for t0 ~= 0
% 4. Slightly increase odeinstep higher than normal if setting delphi = 0

% Known problems:
% ~

% To do:
% 1. Continue solution for booster until landing
% 2. Incorporate wind pressure (redefine FX perhaps)
% 3. Model spheroid earth (eliminate tangent launch plane to earth)
% 4. Energy (need chemical potential energy from propulsion)
% 5. Streamline input process

% Suppress warnings induced by switch options
%#ok<*UNRCH> (unreachable)
%#ok<*DEFNU> (unused func)
%#ok<*NASGU> (unused var)
%#ok<*ALIGN> (misaligned)

% Conversions from [US] to [US] and [US] to [SI]
ft2in = 12; % [ft] --> [in]
fps2mph = 0.681818; % [ft/s] --> [mi/hr]
in2m = 0.0254; % [in] --> [m]
lbf2N = 4.44822; % [lbf] --> [N]
lb2kg = 0.453592; % [lb] --> [kg]
pi2Pa = 6894.76; % [psi] --> [Pa]

% Options
% The vector i:j with 0 < i < j < 9 means show sequences i to j where
% 1 - Rail clearance        3 - Separation Delay      5 - Burnout 2        7 - Parachute
% 2 - Burnout 1             4 - Ignition Delay        6 - Apogee           8 - Touchdown
tableSeqs = 1:6; % i:j - Sequences to display in tables (associated with DispTable)
plotsSeqs = 1:6; % i:j - Sequences to display in plots (associated with Plots)
angleunit = 'rad'; % 'deg'/'rad' - Angle units displayed in plot outputs
USCSIunit = 'SI'; % 'USC'/'SI' - Units in which to display results (US Customary / Metric) (Not yet implemented)
odestats_ = 'on'; % 'on'/'off' - Give ODE solver statistics
oderefine = 8; % Positive integer (affects solution smoothness)
odereltol = 1e-7; % Positive scalar - Relative error tolerance bound
odeabstol = 1e-9; % Positive scalar - Absolute error tolerance bound
odeinstep = 1; % Positive scalar - Initial step size in time
odemxstep = 0.1; % Positive scalar - Maximum permissible step size in time (including inf)
OmitCalcs = false; % true/false - Save time by performing no additional calculations
TinyTable = false; % true/false - Short table showing column names and end of data
FullTable = false; % true/false - Entire table set of all data
DispTable = false; % true/false - Event table
MakePlots = false; % true/false - Show plots
RealTimes = true; % true/false - Show runtime after finishing
ODEsTimes = true; % true/false - Show solver times at each step
HaveStoch = false; % true/false - Consider random (stochastic) forcing
RandStoch = false; % true/false - Consider doubly random forcing

% Stochastic options - Satisifies |n| < 1 and n = 0 implies no stochastics.
StochXmax = 0.05; % Maximum that X (chi) can be. Only associated with HaveStoch.
RandSXmax = 0.05; % Max that n can be sampled to be. Is 0 if RandStoch is false.
% Launch site vector of parameters into solver
geodetLat = 40.9107; % Geodetic latitude of launch site [deg]
localElev = 1191; % Local elevation above sea level [m]
localTemp = 70; % Local surface temperature [F]
towerSpan = 30 * ft2in * in2m; % Length of launch tower [ft] --> [in] --> [m]
% Rocket-specific vector of parameters
allLength = [6.35, 3.71]; % [Overall length during S1, Overall length during S2] [m]
outerDiam = [8.5, 6] * in2m; % [S1 outer diameter, S2 outer diameter] [in] --> [m]
massFuelS = [63.49, 23.58]; % [S1 fuel mass, S2 fuel mass] [kg]
massInitS = [194.43, 65.33]; % [S1 initial mass, S2 initial mass] [kg]
BurnSimFT = {csvread("csvs/S1SL.csv", 1, 0), csvread("csvs/S235.csv", 1, 0)}; % Obtain data from BurnSim csv files
RasAeroCd = {csvread("csvs/H7512.csv", 1, 0), csvread("csvs/H7522.csv", 1, 0)}; % Obtain data from RasAero csv files
specHeata = 1.4; % Freestream ratio of specific heats
chambTemp = [2477, 2477]; % Chamber temperature [K]
specHeatg = [1.25, 1.25]; % [S1 perfect ratio of specific heats, S2 perfect ratio of specific heats] []
diaThroat = [2.5, 1.6] * in2m; % [S1 throat diameter, S2 throat diameter] [in] --> [m]
diaEgress = [6.01, 4.47] * in2m; % [S1 exit diameter, S2 exit diameter] [in] --> [m]
diaFlatDM = [36, 108] * in2m; % [flat drogue diameter, flat main diameter] [in] --> [m]
sepDelays = 2; % Time it takes for second stage to separate from first stage [s]
ignDelays = 0; % Time it takes for second stage to fire after separation [s]
deployAlt = 609; % Altitude of main parachute deployment (2000 ft) [m]
[angleunit, delT, rho0, csvs, tb, tbmx, Ar, At, Ae, Ap, subexitMa, supexitMa, hpc] = CorrectionsAndCalculations();
args.iWant = table(angleunit, OmitCalcs, MakePlots, plotsSeqs, RealTimes, ODEsTimes, HaveStoch, RandStoch);
args.stochasticParams = table(StochXmax, RandSXmax);
args.launchSiteParams = table(geodetLat, localElev, rho0, delT, towerSpan);
args.rocketParams = table(allLength, outerDiam, massFuelS, massInitS, csvs, specHeata, specHeatg, chambTemp, subexitMa, supexitMa, sepDelays, ignDelays, diaFlatDM, deployAlt, tbmx, tb, Ar, At, Ae, Ap, hpc);

% Constants and measured quantities relevant to ODE initial conditions
Lcm = 2.2; % Distance from booster nozzle to center of mass [m]
delr = 3 * in2m + outerDiam(1)/2; % (guess) Perpendicular distance from rail axis to mass center axis [m]
delphi = 5; % Angle made between railing and vertical [deg]
deltheta = 30; % Angle made between x-axis and due East [deg]
veps = eps; % Very small initial velocity [m/s]
args.initParams = table(Lcm, delr, delphi, deltheta, veps);

% ODE initial conditions - Define the state at time t to be xx = [x, y, z, x', y', z']'.
t0 = 0; % No support for nonzero initial time.
tf = 2000; % Final time of propagation [s]
xx0 = [(Lcm * sind(delphi) + delr * cosd(delphi)) * cosd(deltheta);
       (Lcm * sind(delphi) + delr * cosd(delphi)) * sind(deltheta);
       (Lcm * cosd(delphi) - delr * sind(delphi));
       veps * sind(delphi) * cosd(deltheta);
       veps * sind(delphi) * sind(deltheta);
       veps * cosd(delphi)];
if (RealTimes), tic, end % Start timer
% Solve the ODEs with an appropriate solver - for help/pointers, visit the link below
% https://www.mathworks.com/help/matlab/math/troubleshoot-common-ode-problems.html
% Set default options
fsolveoptions = optimoptions('fsolve', 'display', 'none');
options = odeset('reltol', odereltol, 'abstol', odeabstol, 'initialstep', odeinstep, 'maxstep', odemxstep, 'stats', odestats_, 'refine', oderefine, 'events', @odevents);
args.solveropts = table(fsolveoptions, options);
% Create cells to hold information from each sequence
[t, xx, te, xxe] = deal(cell(8, 1));
args.initParams.laste = {zeros(8, 1)};
seqStr(1) = "operating first stage, on launch rail";
seqStr(2) = "operating first stage, rail cleared";
seqStr(3) = "coasting between first and second stages with depleted first stage";
seqStr(4) = "coasting between first and second stages with detached first stage";
seqStr(5) = "operating second stage";
seqStr(6) = "coasting up to apogee";
seqStr(7) = "descending under drogue parachute";
seqStr(8) = "descending under main parachute";
% Begin integrating sequences 1 - 7. Perform seq. 1 separately to initiate recurrence
fprintf("Propagating sequence %1.0f (%s)\n", 1, seqStr(1))
[t{1}, xx{1}, te{1}, xxe{1}, ~] = ode113(@(t, xx, args, seq) odeval(t, xx, args, 1, 1), [t0, tf], xx0, options, args, 1);
stage = 1;
args.initParams.laste{1}(1) = te{1};
seq = 2;
maxseq = max(max(plotsSeqs), max(tableSeqs));
while (seq - 1 ~= maxseq && max(t{seq - 1}) < tf)
    if (seq == 3)
        if (sepDelays ~= 0)
            fprintf("\nPropagating sequence 3 (%s)\n", seqStr(3))
            [t{3}, xx{3}, te{3}, xxe{3}, ~] = ode113(@(t, xx, args, seq) odeval(t, xx, args, 3, stage), [te{2}, tf], xxe{2}, options, args, 3);
        else
            fprintf("\nPropgating sequence 3 - No separation delay, copying final state from sequence 2\n")
            [te{3}, xxe{3}] = deal(te{2}, xxe{2});
        end
        args.initParams.laste{1}(3) = te{3};
        seq = 4;
        stage = 2;
        if (maxseq == 3), break, end
    end
    if (seq == 4)
        if (ignDelays ~= 0)
            fprintf("\nPropagating sequence 4 (%s)\n", seqStr(4))
            [t{4}, xx{4}, te{4}, xxe{4}, ~] = ode113(@(t, xx, args, seq) odeval(t, xx, args, 4, stage), [te{3}, tf], xxe{3}, options, args, 4);
        else
            fprintf("\nPropagating sequence 4 - No ignition delay, copying final state from sequence 3\n")
            [te{4}, xxe{4}] = deal(te{3}, xxe{3});
        end
        args.initParams.laste{1}(4) = te{4};
        seq = 5;
        if (maxseq == 4), break, end
    end
    fprintf("\nPropagating sequence %1.0f (%s)\n", seq, seqStr(seq))
    [t{seq}, xx{seq}, te{seq}, xxe{seq}, ~] = ode113(@(t, xx, args, seq) odeval(t, xx, args, seq, stage), [te{seq - 1}, tf], xxe{seq - 1}, options, args, seq);
    args.initParams.laste{1}(seq) = te{seq};
    seq = seq + 1;
end
% Check if seqs. 4,3 events were used only as holders from seqs. 3,2 and remove if so
if (ignDelays == 0 && maxseq > 3 && te{4} == te{3}), [te{4}, xxe{4}] = deal([], []); end
if (sepDelays == 0 && maxseq > 2 && te{3} == te{2}), [te{3}, xxe{3}] = deal([], []); end

disp("Done propagating. Now calculating.")
% Remove duplicate points from matching [end --> initial] conditions
jj = find(~cellfun(@isempty, t))';
for ii = jj 
    t{ii}(1, :) = [];
    xx{ii}(1, :) = [];
end
% Allocate memory
[x, y, z, vx, vy, vz, rvec, vvec, r, v, rxy, phi, theta, phip, thetap, psi, El, Eld, psid, thetapd, phipd, thetad, phid, ...
    avec, ax, ay, az, a, FT, FD, FG, FP, FN, m, Cd, g, T, p, rho, q, vxv, vyv, vzv, FXx, FXy, FXz, c, Ma, mu, nu, Re] = deal(cell(8, 1));
lenT = cellfun(@length, t); % Size of each sequence
% Calculate (to save time, in parallel if delay ~= 0 by changing for --> parfor)
for ii = jj
    % Check if separation has happened yet
    if (ii > 3)
        stage = 2; % Change stages
        charL = allLength(2); % Use stage 2 characteristic length for Re
    else
        stage = 1;
        charL = allLength(1); % Use stage 1 characteristic length for Re
    end
    % Retrieve state variables from the state (true solutions)
    [x{ii}, y{ii}, z{ii}] = deal(xx{ii}(:, 1), xx{ii}(:, 2), xx{ii}(:, 3));
    [vx{ii}, vy{ii}, vz{ii}] = deal(xx{ii}(:, 4), xx{ii}(:, 5), xx{ii}(:, 6));
    
    % Check if saving time by omitting calculations
    if (OmitCalcs), continue; end
    
    % Relate the state variables back to spherical coordinates in O and O' at each sequence ii
    rvec{ii} = [x{ii}, y{ii}, z{ii}]; % Matrix containing [x, y, z]' at time t [m]
    vvec{ii} = [vx{ii}, vy{ii}, vz{ii}]; % Matrix containing [x', y', z']' at time t [m/s]
    r{ii} = vecnorm(rvec{ii}')'; % Radial distance from origin at time t [m]
    v{ii} = vecnorm(vvec{ii}')'; % Magnitude of velocity at time t [m/s]
    phi{ii} = acos(z{ii} ./ r{ii}); % Tracks how r deviates from +z axis at time t [rad]
    theta{ii} = atan2(y{ii}, x{ii}); % Tracks how r sin(phi) moves in x-y plane at time t [rad]
    phip{ii} = acos(vz{ii} ./ v{ii}); % Tracks how r deviates from +z axis at time t [rad]
    thetap{ii} = atan2(vy{ii}, vx{ii}); % Tracks how r sin(phi) moves in x-y plane at time t [rad]
    psi{ii} = atan2(vecnorm(cross(rvec{ii}', vvec{ii}')), dot(rvec{ii}', vvec{ii}'))'; % Tracks how v deviates from r
    El{ii} = 1.5707963267949 - phi{ii}; % Elevation [rad]
    Eld{ii} = rad2deg(El{ii}); % Elevation [rad] --> [deg]
    psid{ii} = rad2deg(psi{ii}); % psi [rad] --> [deg]
    thetapd{ii} = rad2deg(thetap{ii}); % theta' [rad] --> [deg]
    phipd{ii} = rad2deg(phip{ii}); % phi' [rad] --> [deg]
    thetad{ii} = rad2deg(theta{ii}); % theta [rad] --> [deg]
    phid{ii} = rad2deg(phi{ii}); % phi [rad] --> [deg]
    
    % Additional
    rxy{ii} = sqrt(x{ii}.^2 + y{ii}.^2);
    
    % Calculate atmosphere, body forces, and other parameters experienced during flight
    [FT{ii}, FD{ii}, FG{ii}, FP{ii}, FN{ii}, m{ii}, Cd{ii}, g{ii}, T{ii}, p{ii}, rho{ii}, q{ii}] = deal(zeros(lenT(ii), 1));
    [~, ~, g{ii}, ~, T{ii}, p{ii}, rho{ii}, c{ii}, mu{ii}, nu{ii}, Cd{ii}, m{ii}, FT{ii}, FD{ii}, FG{ii}, FP{ii}, FN{ii}, ~, q{ii}, Ma{ii}]  = odetc(t{ii}, xx{ii}', args, ii, stage);
    
    % Recover flight angles to calculate acceleration
    if (ii == 1) % On rail
        vxv{ii} = sind(delphi) * cosd(deltheta) * ones(lenT(1), 1);
        vyv{ii} = sind(delphi) * sind(deltheta) * ones(lenT(1), 1);
        vzv{ii} = cosd(delphi) * ones(lenT(1), 1);
    else % Free flight
        vxv{ii} = xx{ii}(:, 4) ./ v{ii};
        vyv{ii} = xx{ii}(:, 5) ./ v{ii};
        vzv{ii} = xx{ii}(:, 6) ./ v{ii};
    end
    
    % Precalculate FT - FD for speed
    FTminusFD = FT{ii} - FD{ii};
    
    % Calculate/estimate acceleration and stochastic forces
    if (HaveStoch)
        % Estimate acceleration via finite difference to approximate a derivative
        avec{ii} = diff(vvec{ii}) ./ diff(t{ii}); % Approximate derivative (with 'no' corners) [m/s2]
        avec{ii} = [NaN(1, 3); avec{ii}]; % Insert artificial point to recover full size [m/s2]
        [ax{ii}, ay{ii}, az{ii}] = deal(avec{ii}(:, 1), avec{ii}(:, 2), avec{ii}(:, 3)); % Components [m/s2]
        a{ii} = vecnorm(avec{ii}')'; % Magnitude of estimated acceleration [m/s2]
        % Sequences 1, 5, 6, & 7 have been designated to have no acting stochastics
        if (ii == 1 || ii == 5 || ii == 6 || ii == 7)
            [FXx{ii}, FXy{ii}, FXz{ii}] = deal(zeros(numel(a{ii}), 1));
        else % Estimate stochastic forces
            FXx{ii} = m{ii} .* ax{ii} - (FTminusFD) .* vxv{ii};
            FXy{ii} = m{ii} .* ay{ii} - (FTminusFD) .* vyv{ii};
            FXz{ii} = m{ii} .* az{ii} - (FTminusFD) .* vzv{ii} + FG{ii};
        end
    else
        % Calculate acceleration without finite difference (true full size)
        [FXx{ii}, FXy{ii}, FXz{ii}] = deal(zeros(numel(t{ii}), 1));
        ax{ii} = ((FP{ii} * sind(delphi) - FN{ii} * cosd(delphi)) * cosd(deltheta) + FTminusFD .* vxv{ii} + FXx{ii}) ./ m{ii};
        ay{ii} = ((FP{ii} * sind(delphi) - FN{ii} * cosd(delphi)) * sind(deltheta) + FTminusFD .* vyv{ii} + FXy{ii}) ./ m{ii};
        az{ii} = (FP{ii} * cosd(delphi) + FN{ii} * sind(delphi) + FTminusFD .* vzv{ii} - FG{ii} + FXz{ii}) ./ m{ii};
        avec{ii} = [ax{ii}, ay{ii}, az{ii}];
        a{ii} = vecnorm(avec{ii}')';
    end
    
    % Obtain additional information from the flight using these data
    Re{ii} = v{ii} * charL ./ nu{ii}; % Reynolds number []
end

disp("Done calculating. Now building.")
% Adjust event functions to replace any empty cells with NaN to avoid dynamic table
te(cellfun('isempty', te)) = {NaN};
xxe(cellfun('isempty', xxe)) = {NaN(1, 6)};
% Build a small table to hold the event times and states
evseqs = ["Clearance", "Burnout 1", "Sep. Delay", "Ign. Delay", "Burnout 2", "Apogee", "Parachute", "Landing"]';
YevIDs = table(evseqs(tableSeqs), 'variablenames', "Identifiers");
YevVNs = ["t", "x", "y", "z", "vx", "vy", "vz"];
YevVUs = ["s", "m", "m", "m", "m/s", "m/s", "m/s"];
% Assign titles for minimum and maximum values
YmmIDs = table(["maximum", "minimum"]', 'variablenames', "Identifiers");
% Build the cell (C) table of values
XXtxaC = [t, x, y, z, vx, vy, vz, ax, ay, az]; % Full state - time, pos, vel, acc in x, y, z
XXrvaC = [r, rxy, v, a]; % Magnitudes of pos, vel, acc
FTDGSC = [FT, FD, FG, FP, FN, FXx, FXy, FXz]; % Forces (FX_s very badly estimated)
YvangC = [vxv, vyv, vzv]; % Angles modulating (FT - FD). NOT percentage of v per axis
YangrC = [theta, phi, thetap, phip, psi, El]; % Spherical coordinate angles in radians
YangdC = [thetad, phid, thetapd, phipd, psi, Eld]; % Spherical coordinate angles in degrees
YflytC = [m, Cd, q, Ma, Re]; % Flight parameters
YatmoC = [g, T, p, rho, c, mu, nu]; % Atmospheric parameters
specHeata = [XXtxaC, XXrvaC, FTDGSC, YvangC, YangrC, YangdC, YflytC, YatmoC]; % Create cell table of all values 
% Give names (N) to each column in order of appearance above
XXtxaN = ["t", "x", "y", "z", "vx", "vy", "vz", "ax", "ay", "az"];
XXrvaN = ["r", "rxy", "v", "a"];
FTDGSN = ["FT", "FD", "FG", "FP", "FN", "FXx", "FXy", "FXz"];
YvangN = ["vxv", "vyv", "vzv"];
YangrN = ["theta", "phi", "thetap", "phip", "psi", "El"];
YangdN = ["thetad", "phid", "thetapd", "phipd", "psid", "Eld"];
YflytN = ["m", "Cd", "q", "Ma", "Re"];
YatmoN = ["g", "T", "p", "rho", "c", "mu", "nu"];
YN = [XXtxaN, XXrvaN, FTDGSN, YvangN, YangrN, YangdN, YflytN, YatmoN];
% Give units (U) to each column in order of appearance above
XXtxaU = ["s", "m", "m", "m", "m/s", "m/s", "m/s", "m/s2", "m/s2", "m/s2"];
XXrvaU = ["m", "m", "m/s", "m/s2"];
FTDGSU = ["N", "N", "N", "N", "N", "N", "N", "N"];
YvangU = ["radians", "radians", "radians"];
YangrU = ["radians", "radians", "radians", "radians", "radians", "radians"];
YangdU = ["degrees", "degrees", "degrees", "degrees", "degrees", "degrees"];
YflytU = ["kg", "--", "Pa", "--", "--"];
YatmoU = ["m/s2", "K", "Pa", "kg/m3", "m/s", "kg/ms", "m2/s"];
YU = [XXtxaU, XXrvaU, FTDGSU, YvangU, YangrU, YangdU, YflytU, YatmoU];
% Create tables
if (OmitCalcs), specHeata = specHeata(:, 1:7); YN = YN(:, 1:7); YU = YU(:, 1:7); end % Grab only t and xx if calculations omitted
Ymatrx = cell2mat(specHeata(tableSeqs, :)); % Converted cell to array for use below 
Yplots = cell2table(specHeata(plotsSeqs, :), 'variablenames', YN);
Ysampl = array2table(Ymatrx(end - 5 : end, :), 'variablenames', YN); % Short sample table for reference
Ydsply = array2table(cell2mat(specHeata(tableSeqs, :)), 'variablenames', YN); % Wanted values
YfullY = array2table(cell2mat(specHeata), 'variablenames', YN); % All values
YmnmxY = [YmmIDs, array2table([max(Ymatrx); min(Ymatrx)], 'variablenames', YN)]; % Summary (min/max) values
Yevnts = [YevIDs, array2table(cell2mat([te(tableSeqs), xxe(tableSeqs, :)]), 'variablenames', YevVNs)]; % Event values
% Assign units to tables in order of appearance above
Yplots.Properties.VariableUnits = YU;
Ysampl.Properties.VariableUnits = YU;
Ydsply.Properties.VariableUnits = YU;
YfullY.Properties.VariableUnits = YU;
YmnmxY.Properties.VariableUnits = ["--", YU];
Yevnts.Properties.VariableUnits = ["--", YevVUs];
% Show tables according to inputs - the max/min & event information are always shown last
if (FullTable), disp(YfullY), end
if (DispTable), disp(Ydsply), end
if (TinyTable), disp(Ysampl), end
disp("Maxima and minima of calculated variables")
disp(YmnmxY) % Show
disp("Event states")
disp(Yevnts) % Show

% Add evseqs to args to pass into plotting
args.plottingtools = table(evseqs);

% Remove some unneeded variables from memory to declutter workspace
clear evseqs YevIDs YmmIDs YevVNs YevVUs Ymatrx YN YU
clear XXtxaC XXrvaC FTDGSC YangrC YangdC YflytC YatmoC 
clear XXtxaN XXrvaN FTDGSN YangrN YangdN YflytN YatmoN 
clear XXtxaU XXrvaU FTDGSU YangrU YangdU YflytU YatmoU

% Wrap up
if (MakePlots)
    disp("Done building. Now plotting.")
    plots(Yplots, args)
    disp("Done plotting.")
else
    disp("Done building.")
end

if (~OmitCalcs)
    % Receive a quick estimate for amount of delta-V from launch to full burnout
    [deltaVx, deltaVy, deltaVz, deltaV] = delv([1, 5], t, avec, {FXx, FXy, FXz});
    disp("Launch to full burnout delta-V ~" + deltaV + " m/s.");
end

if (RealTimes), toc, end % End timer

% Reserve space here for possible future plans
disp("Done.")

% Functions
function [angleunit, delT, rho0, csvs, tb, tbmx, Ar, At, Ae, Ap, subexitMa, supexitMa, hpc] = CorrectionsAndCalculations() 
    % Get some unit conversions from base
    pi2Pa = evalin('base', 'pi2Pa');
    lbf2N = evalin('base', 'lbf2N');
    lb2kg = evalin('base', 'lb2kg');
    
    % Get necessary inputs from base
    angleunit = evalin('base', 'angleunit');
    localElev = evalin('base', 'localElev');
    geodetLat = evalin('base', 'geodetLat');
    localTemp = evalin('base', 'localTemp');
    BurnSimFT = evalin('base', 'BurnSimFT');
    RasAeroCd = evalin('base', 'RasAeroCd');
    specHeatg = evalin('base', 'specHeatg');
    outerDiam = evalin('base', 'outerDiam');
    diaThroat = evalin('base', 'diaThroat');
    diaEgress = evalin('base', 'diaEgress');
    diaFlatDM = evalin('base', 'diaFlatDM');
    
    % Correct using " versus '.
    if (~isstring(angleunit)), angleunit = string(angleunit); end
    
    % Obtain difference in local conditions from the standard atmosphere
    [~, ~, ~, ~, T00, ~, ~, ~, ~, ~] = atmos(localElev, geodetLat, 0); % Local surface properties [kg/m3]
    delT = convT(localTemp, 'F', 'K') - T00; % Temperature difference between local and standard at elevation [K]
    [~, ~, ~, ~, ~, ~, rho0, ~, ~, ~] = atmos(localElev, geodetLat, delT); % Local surface density [kg/m3]
    
    % Split data from BurnSimFT and RasAeroCd
    tbpc = {[BurnSimFT{1}(:, 1), BurnSimFT{1}(:, 3)], [BurnSimFT{2}(:, 1), BurnSimFT{2}(:, 3)]}; % Variable chamber pressure in time as csv file [psia]
    tbFT = {[BurnSimFT{1}(:, 1), BurnSimFT{1}(:, 4)], [BurnSimFT{2}(:, 1), BurnSimFT{2}(:, 4)]}; % Variable thrust in time as csv file [lbf]
    tmdt = {[BurnSimFT{1}(:, 1), BurnSimFT{1}(:, 6)], [BurnSimFT{2}(:, 1), BurnSimFT{2}(:, 6)]}; % Variable mass flow rate in time as csv file [lb/s]
    MaCd = {[RasAeroCd{1}(:, 1), RasAeroCd{1}(:, 4:5)], [RasAeroCd{2}(:, 1), RasAeroCd{2}(:, 4:5)]}; % Variable Cd in Mach number as csv file [lbf]
    csvs = {{MaCd{1}, tbpc{1}, tbFT{1}, tmdt{1}}, {MaCd{2}, tbpc{2}, tbFT{2}, tmdt{2}}}; % Combining into one argument | csvs{1} -> S1, csvs{2} -> S2
    
    % Find burn times from csv files
    tbmx = [tbFT{1}(end, 1), tbFT{2}(end, 1)];
    tb = zeros(1, 2); % Burn time vector for stages 1, 2 determined from csv files [s]
    % Calcuable values vector from above section
    Ar = pi * outerDiam.^2 / 4; % [S1 cross-sectional area, S2 cross-sectional area] [m2]
    At = pi * diaThroat.^2 / 4; % [S1 throat area, S2 throat area] [m2]
    Ae = pi * diaEgress.^2 / 4; % [S1 exit area, S2 exit area] [m2]
    Ap = pi * diaFlatDM.^2 / 8; % [Inflated drogue parachute area, inflated main parachute area] [m2]
    fsolveoptions = optimoptions('fsolve', 'display', 'none');
    Y = specHeatg;
    subexitMa = fsolve(@(M) (1 ./ M) .* (2 ./ (Y + 1) .* (1 + (Y - 1) .* M.^2 / 2)).^((1 / 2) * (Y + 1) ./ (Y - 1)) - Ae ./ At, [eps, eps], fsolveoptions); % Solve for exit Mach numbers with initial guesses [0, 0]
    supexitMa = fsolve(@(M) (1 ./ M) .* (2 ./ (Y + 1) .* (1 + (Y - 1) .* M.^2 / 2)).^((1 / 2) * (Y + 1) ./ (Y - 1)) - Ae ./ At, [4, 4], fsolveoptions); % Solve for exit Mach numbers with initial guesses [4, 4]
    clear Y
    hpc{2}{5} = []; % Preallocate necessary positions for interpolating coefficient structures
    for stage = 1:2 
        % Determine tb from BurnSim files (Nonneg. FT at tb = 0)
        tb(stage) = tbFT{stage}(end, 1); % Guess no parasitic data points (chain of FT = 0)
        if (tbFT{stage}(end, 2) == 0) % Find actual burn time otherwise
            [rows, ~] = find(tbFT{stage}(:, 2) == 0); % Obtain rows for which FT = 0
            jj = rows(end); % Start indexing from the end
            while (tbFT{stage}(jj, 2) == 0)
                jj = jj - 1;
            end
            kk = jj;
            if (jj < rows(end)), kk = jj + 1; end
            tb(stage) = tbFT{stage}(kk, 1); % Actual burn time [s]
        end
        
        % Obtain Hermite interpolating polynomial coefficients for chamber pressure, thrust, mass flow, and drag coefficient
        for entry = 1:4
            % Account for unit conversions
            switch entry, case 1, unitconv = 1; case 2, unitconv = pi2Pa; case 3, unitconv = lbf2N; case 4, unitconv = -lb2kg; end
            switch entry
                case 1
                    hpc{stage}{entry} = pchip(csvs{stage}{entry}(:, 1), [csvs{stage}{entry}(:, 2), csvs{stage}{entry}(:, 3)]');
                otherwise
                    hpc{stage}{entry} = pchip(csvs{stage}{entry}(:, 1), csvs{stage}{entry}(:, 2) * unitconv);
            end
        end
        hpc{stage}{entry + 1} = pchip(tmdt{stage}(:, 1), cumtrapz(tmdt{stage}(:, 1), tmdt{stage}(:, 2)) * (-lb2kg));
        clear unitconv
    end % Set up fixed numbers
end
function [Rlat, glat, g, Zg, T, p, rho, c, mu, nu] = atmos(z, lat, delT0) 
    % Constants - Note lat is the geodetic latitude in degrees
    Req = 6378137.14; % Earth's equatorial radius [m]
    Rpo = 6356752.22; % Earth's polar radius [m]
    % Calculate the gravitational constant at geodetic latitude lat
    glat = 9.780356 * (1 + 0.0052885 * sind(lat)^2 - 0.0000059 * sind(2 * lat).^2);
    % Calculate Earth's WGS84 radius at geodetic latitude lat
    Rlat = sqrt(((Req^2 * cosd(lat)).^2 + (Rpo^2 * sind(lat)).^2) ./ ((Req * cosd(lat)).^2 + (Rpo * sind(lat)).^2));
    % Calculate local gravitational acceleration at geodetic latitude lat and altitude z above WGS84 ellipsoid
    g = glat .* (Rlat ./ (Rlat + z)).^2;
    % Convert altitude (geometric) z to altitude (geopotential) ZG [km]
    Zg = (Rlat .* z ./ (Rlat + z)) / 1000;
    % Acquire length of z to prepare in the case if z is nonscalar (i.e. vector-valued)
    lenz = numel(z);
    [T, p] = deal(NaN(lenz, 1)); % Begin outputs (T, p) as NaNs
    % http://www.braeunig.us/space/atmmodel.htm (2014) (WMO 1976 US Standard)
    for ii = 1:lenz % Iteratate since z is possibly a vector
        Zgii = Zg(ii); % Initialize 
        if (Zgii <= 11)
            Tii = 288.15 - 6.5 * Zgii;
            pii = 101325 * (288.15 / Tii)^(-34.1632 / 6.5);
        elseif (Zgii <= 20)
            Tii = 216.65;
            pii = 22632.06 * exp(-34.1632 * (Zgii - 11) / Tii);
        elseif (Zgii <= 32)
            Tii = 196.65 + Zgii;
            pii = 5474.889 * (216.65 / (216.65 + (Zgii - 20)))^(34.1632);
        elseif (Zgii <= 47)
            Tii = 139.05 + 2.8 * Zgii;
            pii = 868.0187 * (228.65 / (228.65 + 2.8 * (Zgii - 32)))^(34.1632 / 2.8);
        elseif (Zgii <= 51)
            Tii = 270.65;
            pii = 110.9063 * exp(-34.1632 * (Zgii - 47) / Tii);
        elseif (Zgii <= 71)
            Tii = 413.45 - 2.8 * Zgii;
            pii = 66.93887 * (270.65 / (270.65 - 2.8 * (Zgii - 51)))^(34.1632 / -2.8);
        elseif (Zgii <= 84.852)
            Tii = 356.65 - 2 * Zgii;
            pii = 3.95642 * (214.65 / (214.65 - 2 * (Zgii - 71)))^(34.1632 / -2);
        elseif (Zgii > 84.852) % continue modeling for very high geocentric altitudes z < 1000 km.
            ziikm = z(ii) / 1000; % Altitude [m] --> [km]
            ziikmv = [ziikm^4, ziikm^3, ziikm^2, ziikm, 1];
            if (ziikm <= 91)
                P = [0.000000, 	2.159582E-06, 	-4.836957E-04, 	-0.1425192, 	13.47530]';
                R = [0.000000, 	-3.322622E-06, 	9.111460E-04, 	-0.2609971, 	5.944694]';
                Tii = 186.8673;
            elseif (ziikm <= 100)
                P = [0.000000, 	3.304895E-05, 	-0.009062730, 	0.6516698, 	-11.03037]';
                R = [0.000000, 	2.873405E-05, 	-0.008492037, 	0.6541179, 	-23.62010]';
                Tii = 263.1905 - 76.3232 * sqrt(1 - ((ziikm - 91) / -19.9429)^2);
            elseif (ziikm <= 110)
                P = [0.000000, 	6.693926E-05, 	-0.01945388, 	1.719080, 	-47.75030]';
                R = [-1.240774E-05, 	0.005162063, 	-0.8048342, 	55.55996, 	-1443.338]';
                Tii = 263.1905 - 76.3232 * sqrt(1 - ((ziikm - 91) / -19.9429)^2);
            elseif (ziikm <= 120)
                P = [0.000000, 	-6.539316E-05, 	0.02485568, 	-3.223620, 	135.9355]';
                R = [0.00000, 	-8.854164E-05, 	0.03373254, 	-4.390837, 	176.5294]';
                Tii = 240 + 12 * (ziikm - 110);
            elseif (ziikm <= 150)
                P = [2.283506E-07, 	-1.343221E-04, 	0.02999016, 	-3.055446, 	113.5764]';
                R = [3.661771E-07, 	-2.154344E-04, 	0.04809214, 	-4.884744, 	172.3597]';
            elseif (ziikm <= 200)
                P = [1.209434E-08, 	-9.692458E-06, 	0.003002041, 	-0.4523015, 	19.19151]';
                R = [1.906032E-08,	-1.527799E-05, 	0.004724294, 	-0.6992340, 	20.50921]';
            elseif (ziikm <= 300)
                P = [8.113942E-10, 	-9.822568E-07, 	4.687616E-04, 	-0.1231710, 	3.067409]';
                R = [1.199282E-09, 	-1.451051E-06, 	6.910474E-04, 	-0.1736220, 	-5.321644]';
            elseif (ziikm <= 500)
                P = [9.814674E-11, 	-1.654439E-07, 	1.148115E-04, 	-0.05431334, 	-2.011365]';
                R = [1.140564E-10, 	-2.130756E-07, 	1.570762E-04, 	-0.07029296, 	-12.89844]';
            elseif (ziikm <= 750)
                P = [-7.835161E-11, 	1.964589E-07, 	-1.657213E-04, 	0.04305869,      -14.77132]';
                R = [8.105631E-12, 	-2.358417E-09, 	-2.635110E-06, 	-0.01562608, 	-20.02246]';
            elseif (ziikm <= 1000)
                P = [2.813255E-11, 	-1.120689E-07, 	1.695568E-04, 	-0.1188941,      14.56718]';
                R = [-3.701195E-12, 	-8.608611E-09, 	5.118829E-05, 	-0.06600998, 	-6.137674]';
            else
                disp("No atmospheric data for altitudes exceeding 1000 km.")
                P = nan(5, 1);
                R = nan(5, 1);
            end
            if (~exist('Tii', 'var')) % z is between 120 km and 1000 km.
                xi = (ziikm - 120) * (6356.766 + 120) / (6356.766 + ziikm);
                Tii = 1000 - 640 * exp(-0.01875 * xi);
            end
            pii = exp(ziikmv * P);
            rhoii = exp(ziikmv * R);
        end
        % Add intermediate results of T and P to final result
        T(ii) = Tii;
        p(ii) = pii;
        clear Tii
    end
    % Account for temperature corrections (pressure stays constant)
    T = T + delT0; % Corrected temperature
    rho = p ./ (287.053 * T); % Corrected density

    % General atmospheric results independent of model
    c = sqrt(401.87 * T); % Local freestream speed of sound
    mu = 0.006581718 * (T / 273.15).^(3 / 2) ./ (T + 110.4); % Dynamic viscosity [kg/(m s)] (6.581718e-3 = C1)
    nu = mu ./ rho; % Kinematic viscosity [m2/s]
end  % Atmospheric model
function Tconv = convT(T, unitT, unitTconv) 
    switch unitT
        case 'F'
            switch unitTconv
                case 'F'
                    Tconv = T;
                case 'C'
                    Tconv = (5 / 9) * (T - 32);
                case 'K'
                    Tconv = (5 / 9) * (T - 32) + 273.15; % (convT(convT(T, 'F', 'C'), 'C', 'K'))
                case 'R'
                    Tconv = T + 459.67;
                otherwise
                    disp("Invalid Temperature return. Please convert to 'F', 'C', 'K', or 'R'.")
            end
        case 'C'
            switch unitTconv
                case 'F'
                    Tconv = (9 / 5) * T + 32;
                case 'C'
                    Tconv = T;
                case 'K'
                    Tconv = T + 273.15;
                case 'R'
                    Tconv = (9 / 5) * T + 491.67; % (convT(convT(T, 'C', 'F'), 'F', 'R');
                otherwise
                    disp("Invalid Temperature return. Please convert to 'F', 'C', 'K', or 'R'.")
            end
        case 'K'
            switch unitTconv
                case 'F'
                    Tconv = (9 / 5) * T - 459.67; % (convT(convT(T, 'K', 'C'), 'C', 'F'))
                case 'C'
                    Tconv = T - 273.15;
                case 'K'
                    Tconv = T;
                case 'R'
                    Tconv = (9 / 5) * T; % (convT(convT(T, "K", "C"), "C", "R"))
                otherwise
                    disp("Invalid Temperature return. Please convert to 'F', 'C', 'K', or 'R'.")
            end
        case 'R'
            switch unitTconv
                case 'F'
                    Tconv = T - 459.67;
                case 'C'
                    Tconv = (5 / 9) * T - 273.15; % (convT(convT(T, 'R', 'F'), 'F', 'C'))
                case 'K'
                    Tconv = (5 / 9) * T; % (convT(convT(T, 'R', 'C'), 'C', 'K'))
                case 'R'
                    Tconv = T;
                otherwise
                    disp("Invalid Temperature return. Please convert to 'F', 'C', 'K', or 'R'.")
            end
        otherwise
            disp("Invalid Temperature input. Please convert from 'F', 'C', 'K', or 'R'.")
    end
end  % Temperature conversion
function pexit = nozzle(t, args, stage, pinf) 
    % Calorically perfect gas with choked flow but no assumptions about a shock in the nozzle
    hpc = args.rocketParams.hpc;
    pc = ppval(hpc{stage}{2}, t); % Find values for chamber pressure at *>burn<* time t
    Y = args.rocketParams.specHeatg(stage); % Ratio of specific heats for a calorically perfect gas []
    subexitMa = args.rocketParams.subexitMa(stage); % Subsonic Mach number at exit
    supexitMa = args.rocketParams.supexitMa(stage); % Supersonic Mach number at exit (no shocks)
    Msmooth = [subexitMa, supexitMa]; % Sub- and supersonic Mach numbers at exit in presence of no waves
    numelt = numel(t);
    opts = args.solveropts.fsolveoptions;
    At = args.rocketParams.At(stage); % Throat area [m2]
    Ae = args.rocketParams.Ae(stage); % Exit area [m2]
    N = 10; % Number of area points to consider
    
    % Assume first that the exit pressure is the back pressure (atmospheric presure) and Ma = 0 (no flow)
    pexit = pinf;
    Maexit = zeros(numelt, 1);
    
    % Find pressure ratios (back / chamber) to determine if there is a shock in the nozzle (pratio)
    % 1 - Choked subsonic          3 - Design
    % 2 - Exit Normal shock        4 - Actual
    pratio{1} = (1 + (Y - 1) * Msmooth(1)^2 / 2)^(Y / (1 - Y)) * ones(numelt, 1); % Subsonic pressure ratio (pb / p0)|sub (choked)
    pratio{3} = (1 + (Y - 1) * Msmooth(2)^2 / 2)^(Y / (1 - Y)) * ones(numelt, 1); % Supersonic pressure ratio (pb / p0)|sup (design)
    pratio{2} = (1 + 2 * Y / (Y + 1) * (Msmooth(2)^2 - 1)) * pratio{3}; % Pressure ratio for shock at exit (pb / p0)|shk (exit shock)
    pratio{4} = pinf ./ pc; % Current pressure ratio (pb / p0)|now
    % Find the flow regime throughout the nozzle at time t (ind)
    % 1 - Subsonic flow        3 - Oblique shock (Overexpanded)
    % 2 - Normal shock         4 - Expansion fan (Underexpanded)
    [ind{1}, ~] = find(pratio{1} <= pratio{4}); % Find indices for which flow is fully subsonic (NOT choked)
    [ind{2}, ~] = find(pratio{2} < pratio{4} & pratio{4} <= pratio{1}); % Indices for which there is shock in the nozzle
    [ind{3}, ~] = find(pratio{3} < pratio{4} & pratio{4} < pratio{2}); % Indices for which flow is fully supersonic in the nozzle but with oblique shock separating pe from pb
    [ind{4}, ~] = find(pratio{4} <= pratio{3}); % Indices for which flow is fully supersonic in the nozzle but with expansion fans separating pe from pb
    
    % Find which, if any, indices are not empty and which are
    nn = ~cellfun(@isempty, ind); % Indices that are not empty (there is a flow regime at time t here)
    jj = find(nn); % Positions in kk that are 1 (i.e. not empty)
    [xMa, xTTcRatio, xppcRatio, xrhorhocRatio] = deal(NaN(numelt, N));
    [Xshock, beta, theta, thetaTA] = deal(NaN(numelt, 1));
    for ii = jj
        kk = ind{ii}; % Positions in fill results into for current ii
        numelkk = numel(kk);
        switch ii % Switch 1 - 4 based on definitions for ind (subsonic flow - expansion fans)
            case 1 % Subsonic flow
                pexit(kk) = pinf(kk); % Exit pressure [Pa]
                Maexit(kk) = 2 / (Y - 1) * (pratio{4}(kk).^(1/Y - 1) - 1); % Exit Mach number []
                if (numelt > 1)
                    % Calculate required A* for this flow to be choked
                    Astar = Maexit(kk) .* (2 / (Y + 1) * (1 + (Y - 1) * Maexit(kk).^2 / 2)).^(-(1 / 2) * (Y + 1) / (Y - 1)) * Ae; % A* - Area at which flow would be sonic [m2]
                    Mathroat = fsolve(@(M)  (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - At ./ Astar, 0.2 * ones(numelkk, 1), opts); % Throat mach number []
                    % Determine ratios along the nozzle (with diverging section monotonically increasing) at stations x
                    A = linspace(At, Ae, N); % Areas to consider in evaluating isentropic ratios - Assumes monotonically increasing nozzle area
                    X = (A - At) / (Ae - At); % Location along nozzle for which ratios are evaluated (0 < x < 1)
                    for ll = kk'
                        xMa(ll, :) = fsolve(@(M) (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - A / Astar(ll), 0.3*ones(1, N), opts); % Mach number at station
                    end
                    xTTcRatio = (1 + (Y - 1) * xMa(kk, :).^2 / 2).^(-1); % Temperature ratio at stations
                    xppcRatio = xTTcRatio(kk, :).^(Y / (Y - 1)); % Pressure ratio at stations
                    xrhorhocRatio = xppcRatio(kk, :).^(1 / Y); % Density ratio at stations
                    % Solutions: xMa, xTTcRatio, xppcRatio, xrhorhocRatio
                end
            case 2 % Supersonic flow terminated by normal shock to subsonic flow in nozzle
                pexit(kk) = pinf(kk);
                Maexit(kk) = sqrt(sqrt(1 / (Y - 1)^2 + 2 / (Y - 1) * (2 / (Y + 1))^((Y + 1) / (Y - 1)) * (((At / Ae) ./ pratio{4}(kk))).^2) - 1 / (Y - 1)); % Exit Mach number with shock in nozzle []
                if (numelt > 1) % Find shock position for data analysis (not relevant for obtaining the state solution)
                    % Form pressure ratios
                    p0eperatio = (1 + (Y - 1) * Maexit(kk).^2 / 2).^(Y / (Y - 1)); % Exterior (outside of nozzle) pressure ratio
                    p0epefunc = @(M) (1 + 2 * Y / (Y + 1) *(M.^2 - 1)).^(1 + Y / (1 - Y)) .* ((2 + (Y - 1) * M.^2) ./ ((Y + 1) * M.^2)).^(Y / (1 - Y)); % Relate upstream Mach no. to exterior pratio
                    p0epcfunc = @(M) p0epefunc(M) .* pc(kk) ./ pinf(kk); % Exit to chamber stagnation pressure ratios []
                    % Solve for Mach number just upstream of the normal shock
                    MaUSshock = fsolve(@(M) p0epcfunc(M) - p0eperatio, 5 * ones(numelkk, 1), opts); % []
                    % Find the area of shock, fully determining its location in the nozzle
                    Ashock = (2 / (Y + 1) * (1 + (Y - 1) * MaUSshock.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ MaUSshock * At; % [m2]
                    % Express the shock location in the div. section as nondimensionalized position x such that d(t, e) = 1
                    Xshock(kk) = (Ashock - At) / (Ae - At); % Shock location (nondimensionalized between 0 and 1) []
                    % Discretize other areas throughout the diverging portion of the nozzle
                    A = linspace(At, Ae, N);
                    % Express in terms of position
                    X = (A - At) / (Ae - At); % Valid since A = A(x) (i.e. X == linspace(0, 1, N) is true)
                    % Find indices for locations in the nozzle
                    [XUS, tUS] = find(X' < Xshock(kk)');
                    [XDS, tDS] = find(Xshock(kk)' <= X');
                    % Allocate storage to hold storage of solutions at each t
                    [xMaUS, xMaDS] = deal(cell(numelkk, 1));
                    % Solve for supersonic Ma number once since flow is choked
                    sol = fsolve(@(M) (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - A / At, mean(MaUSshock)*ones(1, N), opts);
                    % Mach number solutions upstream of shock (M > 1)
                    for ll = unique(tUS)' % Solve equations at each t (allows to combine with sub. & is generally faster than solving once for all t)
                        [xind, ~] = find(~(tUS - ll)); % Obtain relevant indices for which this step is valid
                        tXUS = XUS(xind); % Find those XDS values corresponding to ll in tDS (i.e. 1xn vector for positions)
                        % Prescribe supersonic solution to correct indices - Mach numbers at stations upstream of shock []
                        xMaUS{ll} = sol(tXUS);
                    end
                    % Find A* in subsonic section required for sonic flow to occur
                    AstarDS = Maexit(kk) .* (2 / (Y + 1) * (1 + (Y - 1) * Maexit(kk).^2 / 2)).^(-(1 / 2) * (Y + 1) / (Y - 1)) * Ae; % A* downstream of shock [m2]
                    % Mach number solutions downstream of shock (M < 1)
                    for ll = unique(tDS)' % Solve equations at each t since AstarDS is generally changing in time as the shock moves
                        [xind, ~] = find(~(tDS - ll)); % Obtain relevant indices for which this step is valid
                        tXDS = XDS(xind); % Find those XDS values corresponding to ll in tDS (i.e. 1xn vector for positions)
                        % Solve - Mach numbers at stations downstream of shock []
                        xMaDS{ll} = fsolve(@(M) (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - A(tXDS) / AstarDS(ll)', mean(Maexit(kk))*ones(size(tXDS))', opts);
                    end
                    xMacell = [xMaUS, xMaDS]; % Concatenate together
                    for mm = kk'
                        xMa(mm, 1:N) = [xMacell{mm - min(kk) + 1, 1}, xMacell{mm - min(kk) + 1, 2}]; % Populate the matrix
                    end
                    % Form time-varying ratios
                    xTTcRatio(kk, :) = 1 ./ (1 + (Y - 1) * xMa(kk, :).^2 / 2); % Temperature ratio at stations
                    xppcRatio(kk, :) = xTTcRatio(kk, :).^(Y / (Y - 1)); % Pressure ratio at stations
                    xrhorhocRatio(kk, :) = xppcRatio(kk, :).^(1 / Y); % Density ratio at stations
                    % Solutions: Xshock, xMa, xTTcRatio, xppcRatio, xrhorhocRatio
                end
            case 3 % Supersonic flow with oblique shock eminating from nozzle exit
                pexit(kk) = pratio{3}(kk) .* pc(kk); % Exit pressure is the same as design
                Maexit(kk) = Msmooth(2); % Supersonic solution from area-Mach relation
                if (numelt > 1)
                    % Approximate shock angle as a linear function of the back chamber pressure ratio
                    perc = (pratio{4}(kk) - pratio{3}(kk)) ./ (pratio{2}(kk) - pratio{3}(kk));
                    % Oblique shock angle relative to exit plane, where along the plane is 0 and normal to it is 90 deg
                    beta(kk) = pi / 2 - (pi / 4) * (1 - perc); % Linearly approximated oblique shock angle relative to exit plane [rad]
                    % Flow turn angle - will be negative as the turn is concave. No quadrant check necessary
                    theta(kk) = atan(2 * cot(beta(kk)) .* (Maexit(kk) .* sin(beta(kk)).^2 - 1) ./ (Maexit(kk).^2 .* (Y + cos(2 * beta(kk))) + 2)); % Flow turn angle [rad]
                    % Determine ratios along the nozzle (with diverging section monotonically increasing) at stations x
                    A = linspace(At, Ae, N); % Areas to consider in evaluating isentropic ratios - Assumes monotonically increasing nozzle area
                    X = (A - At) / (Ae - At); % Location along nozzle for which ratios are evaluated (0 < X < 1)
                    xMavec = fsolve(@(M) (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - A ./ At, Msmooth(2)*ones(1, N), opts); % Mach number at station
                    xMa(kk, :) = ones(numelkk, 1) * xMavec;
                    xTTcRatio(kk, :) = (1 + (Y - 1) * xMa(kk, :).^2 / 2).^(-1); % Temperature ratio at stations
                    xppcRatio(kk, :) = xTTcRatio(kk, :).^(Y / (Y - 1)); % Pressure ratio at stations
                    xrhorhocRatio(kk, :) = xppcRatio(kk, :).^(1 / Y); % Density ratio at stations
                    % Solutions: beta, theta, xMa, xTTcRatio, xppcRatio, xrhorhocRatio
                end
            case 4 % Supersonic flow with expansion fan eminating from nozzle exit
                pexit(kk) = pratio{3}(kk) .* pc(kk); % Exit pressure is the same as design
                Maexit(kk) = Msmooth(2); % Supersonic solution from area-Mach relation
                if (numelt > 1)
                    % Find turning angle - expansion fans are isentropic so total pressure is const.
                    PMaus = sqrt((Y + 1) / (Y - 1)) * atan(sqrt((Y - 1) / (Y + 1) * (Maexit(kk).^2 - 1))) - atan(sqrt(Maexit(kk).^2 - 1)); % Prandtl-Meyer angle immediately before expansion fans
                    Mads = sqrt((2 / (Y - 1)) * ((pinf(kk) ./ pc(kk)).^(1/Y - 1) - 1)); % Mach number through fan []
                    PMads = sqrt((Y + 1) / (Y - 1)) * atan(sqrt((Y - 1) / (Y + 1) * (Mads.^2 - 1))) - atan(sqrt(Mads.^2 - 1)); % Prandyl-Meyer angle immediately after expansion fans
                    thetaTA(kk) = PMads - PMaus; % Turning angle of flow from the nozzle exit through the fan [rad]
                    % Determine ratios along the nozzle (with diverging section monotonically increasing) at stations X
                    A = linspace(At, Ae, N); % Areas to consider in evaluating isentropic ratios - Assumes monotonically increasing nozzle area
                    X = (A - At) / (Ae - At); % Location along nozzle for which ratios are evaluated (0 < X < 1)
                    xMavec = fsolve(@(M) (2 / (Y + 1) * (1 + (Y - 1) * M.^2 / 2)).^((1 / 2) * (Y + 1) / (Y - 1)) ./ M - A ./ At, Msmooth(2)*ones(1, N), opts); % Mach number at station
                    xMa(kk, :) = ones(numelkk, 1) * xMavec;
                    xTTcRatio(kk, :) = (1 + (Y - 1) * xMa(kk, :).^2 / 2).^(-1); % Temperature ratio at stations
                    xppcRatio(kk, :) = xTTcRatio(kk, :).^(Y / (Y - 1)); % Pressure ratio at stations
                    xrhorhocRatio(kk, :) = xppcRatio(kk, :).^(1 / Y); % Density ratio at stations
                    % Solutions: thetaTA, xMa, xTTcRatio, xppcRatio, xrhorhocRatio
                end
        end
    end
end  % Nozzle characteristics
function Cd = coefficients(pinf, Tinf, Minf, RoughnessCoeff, args) 
    % Implementation of J. Barrowman dissertation (1967)
    
    % Work in progress: Not used live in code yet
    
    % To Do for future fixes:
        % Dynamically update N until tangentConeX is automatically sorted (no overlaps between tangent points and regions)
        % Check if Mach number on front cone is subsonic. If so, then no P-M fans
    Rs = RoughnessCoeff;
    Y = args.rocketParams.specHeata;
    % Begin with second-order shock expansion theory on nose cone
    noseConeR = 3 * 0.0254; % Exposed nose cone base radius [in] --> [m]
    noseConeL = 30 * 0.0254; % Exposed nose cone length [in] --> [m]
    noseConeN = 5; % N tangent points (N - 1 section - 1 circular cone and N - 2 conical frustums)
    noseConeX = linspace(eps, noseConeL - eps, noseConeN); % Discretized x locations along cone to be used as tangent points for frustums [m]
    % von Karman ogive nose cone (specific to H75)
    vnKarmanT = @(X) acos(1 - 2 * X / noseConeL); % Parameter used in defining von Karman ogive nose cones []
    HaackC = 0; % Haack's C value determining various kinds of vK ogives []
    noseConer = @(X) noseConeR * sqrt((1 / pi) * (vnKarmanT(X) - sin(2 * vnKarmanT(X)) / 2 + HaackC * sin(vnKarmanT(X)).^3)); % Nose cone radius [m]
	noseCones = @(X) (noseConeR ./ sqrt((noseConeL - X) .* X)) .* (2 + 3 * HaackC * cos(vnKarmanT(X))) .* sin(vnKarmanT(X)).^2 ./ sqrt(2 * pi * (2 * vnKarmanT(X) - sin(2 * vnKarmanT(X)) + 2 * HaackC * sin(vnKarmanT(X)).^3)); % Nose cone slope []
%     % Obtain nose cone radius and slope at each station X
%     noseConeR = noseConer(noseConeX); % Radius from centerline [m]
%     noseConeS = noseCones(noseConeX); % Slope [m/m]
%     % Apply a correction so that all elements, particularly the first, of noseConeS are finite
%     nanX = isnan(noseConeS);
%     if (~all(nanX)), noseConeS(nanX) = spline(noseConeX(~nanX), noseConeS(~nanX), noseConeX(nanX)); end
%     % Allocate memory
%     [tangentConeX, tangentConeT, tangentConeD] = deal(NaN(1, noseConeN - 1));
%     % Find positions on the longitudinal axis where tangent lines intersect (xi)
%     for ii = 1 : noseConeN - 1
%         % Note: (Already in To Do list) Very little protection offered from cutting inside the actual nose cone (noseConeX(ii) < noseConeX(ii + 1) < tangentConeX(ii))
%         % X-coordinate of adjacent tangent line intersection
%         tangentConeX(ii) = - (noseConeR(ii + 1) - noseConeS(ii + 1) * noseConeX(ii + 1) - (noseConeR(ii) - noseConeS(ii) * noseConeX(ii))) / (noseConeS(ii + 1) - noseConeS(ii));
%         
%         % Equivalent cone semi-vertex angle (T, theta) and turning angle (D, delta)
%         if (ii == 1) % Initial circular cone
%             tangentConeT(1) = noseConeS(1); % Equivalent cone semi-vertex angle [rad]
%             tangentConeD(1) = -noseConeS(1); % Surface turning angle (P-M angle) [rad]
%         elseif (ii == 2) % First circular frustum after circular cone
%             tangentConeT(2) = noseConeS(2); % Equivalent cone semi-vertex angle [rad]
%             tangentConeD(2) = tangentConeT(1) - tangentConeT(2); % Surface turning angle (P-M angle) [rad]
%         else % All frustums wedged between other frustums
%             tangentConeT(ii) = noseConeS(ii); % Equivalent cone semi-vertex angle [rad]
%             tangentConeD(ii) = tangentConeT(ii - 1) - tangentConeT(ii); % Surface turning angle (P-M angle) [rad]
%         end
%     end
    
    options = optimoptions('fsolve', 'TolFun', 1e-10, 'display', 'none');
    % Approximate the nose cone with a series of tangent cones - begin algorithm in the middle at x = L/2
    xTangent(1) = noseConeL / 2; % x-location of first tangent point
    yTangent(1) = noseConer(xTangent(1)); % y-location of first tangent point
    mTangent(1) = noseCones(xTangent(1)); % slope at first tangent point

    % (x/y)Connect points are those associated with corners in tangent approximation
    xConnect = []; % x-locations where tangent lines intersect
    yConnect = []; % y-locations where tangent lines intersect
    
    % Options
    maxm = 50; % Maximum permissible tangent slope []
    tol = 0.5e-1; % Tolerance for fitting tangent lines to curve
        
    % Curve is concave so the max slope occurs at x = 0
    maxSlope = noseCones(0);
    if (isnan(maxSlope) || isinf(maxSlope)), maxSlope = noseCones(eps * noseConeL); end

    % Set flags to begin and end the algorithm
    noseReached = false;
    baseReached = false;
    % Set additional flag to check for conditions that clearly end the algorithm
    flag = false;
    % Going towards the nose
    while (~noseReached)
        % Define points and function to be evaluated when going leftwards from L/2
        xSamplePoints = linspace(xTangent(1), 0, 1000);
        numelx = numel(xSamplePoints);
        tangentr = @(x) mTangent(1) * (x - xTangent(1)) + yTangent(1);
        divfunc = @(x) [tangentr(x)./ noseConer(x) - 1, atan(tangentr(x) / x) - atan(noseConer(x) / x)];
        for ii = 1:numelx
            div = divfunc(xSamplePoints(ii)); % Divergence away from actual nose cone in terms of measured area
            if (isnan(div(2))), div(2) = atan(tangentr(eps^2) / eps^2) - atan(noseConer(eps^2) / eps^2); end
            if (atand(mTangent(1)) > 50 || (ii == numelx && norm(div) < tol)) % Reached the end but still within the allowance
                % Connect from current tangent point to tip with secant line
                xConnect = [0, xConnect];
                yConnect = [0, yConnect];
                noseReached = true; % Exit condition
                break
            elseif (norm(div) >= tol) % Exceeded the allowance, create a new point and find new tangent point
                % Create point here that next tangent line has to connect with
                NewxConnect = fsolve(@(x) norm(divfunc(x)) - tol, xSamplePoints(ii), options); % Solution guaranteed
                NewyConnect = mTangent(1) * (NewxConnect - xTangent(1)) + yTangent(1);
                xConnect = [NewxConnect, xConnect];
                yConnect = [NewyConnect, yConnect];
                % Solve for next xTangent
                NewxTangent = fsolve(@(X) noseCones(X) .* (NewxConnect - X) + noseConer(X) - NewyConnect, eps * noseConeL, options);
                if (NewxTangent > xConnect(1))
                    error("New tangent point cannot be behind connection point")
                end
                if (NewxTangent < 0)
                    flag = true;
                    break
                end
                xTangent = [NewxTangent, xTangent];
                yTangent = [noseConer(NewxTangent), yTangent];
                mTangent = [noseCones(NewxTangent), mTangent];
                break
            end
        end
    end
    flag = false;
    while (~baseReached)
        % Define points and function to be evaluated when going leftwards from L/2
        xSamplePoints = linspace(xTangent(end), noseConeL, 1000);
        numelx = numel(xSamplePoints);
        tangentr = @(x) mTangent(end) * (x - xTangent(end)) + yTangent(end);
        divfunc = @(x) [tangentr(x)./ noseConer(x) - 1, tangentr(x) - noseConeR];
        for ii = 1:numelx
            div = divfunc(xSamplePoints(ii)); % Divergence away from actual nose cone in terms of measured area
            if ((ii == numelx && norm(div) < tol) || flag) % Reached the end but still within the allowance
                % Connect from current tangent point to tip with secant line
                xConnect = [xConnect, noseConeL];
                yConnect = [yConnect, noseConeR];
                baseReached = true; % Exit condition
                break
            elseif (norm(div) >= tol) % Exceeded the allowance, create a new point and find new tangent point
                % Create point here that next tangent line has to connect with
                NewxConnect = fsolve(@(x) norm(divfunc(x)) - tol, xSamplePoints(ii), options); % Solution guaranteed
                NewyConnect = mTangent(end) * (NewxConnect - xTangent(end)) + yTangent(end);
                xConnect = [xConnect, NewxConnect];
                yConnect = [yConnect, NewyConnect];
                % Solve for next xTangent
                NewxTangent = fsolve(@(x) noseCones(x) .* (NewxConnect - x) + noseConer(x) - NewyConnect, (NewxConnect + noseConeL) / 2, options);
                if (NewxTangent < xConnect(end))
                    error("New tangent point cannot be ahead of connection point. Try adjusting tol.")
                end
                if (abs(NewxTangent) > noseConeL)
                    flag = true;
                    break
                end
                xTangent = [xTangent, NewxTangent];
                yTangent = [yTangent, noseConer(NewxTangent)];
                mTangent = [mTangent, noseCones(NewxTangent)];
                break
            end
        end
    end

    numCones = numel(xConnect) - 1;
    % Apply correction to mTangent at the beginning and end
    mTangent(1) = yConnect(2) / xConnect(2);
    mTangent(end) = (noseConeR - yTangent(end - 1)) / (noseConeL - xTangent(end - 1));

%     plot(xConnect, yConnect)
%     hold on
%     axis equal
%     r = noseConer(linspace(eps, noseConeL - eps, 100e4));
%     plot(linspace(eps, noseConeL - eps, 100e4), r)
%     grid on
%     grid minor
    
    % Obtain pressure fields over tangent cone surfaces
    if (Minf > 1)
        pc = zeros(1, numCones); % Constant pressure on equivalent cones tangent to the body
        results = zeros(3, numCones);
        S = load('mats/TheTaylorMaccollGrid500x500FilledAndExtended.mat'); % TODO: PULL THIS LINE OUT AND PUT IT IN INPUTS SECTION (RUN ONCE)
        TayMacGrid = S.TheTaylorMaccollGrid;
        % Element (1, 1, :) of TayMacGrid is NaN by design - Split into separate matrices for less confusion in indexing
        TayMacOc = TayMacGrid(1, 2 : end, 1); % Cone angles used in the grid (used for column indices only). dim(TayMacOc) = 1 x (size(TayMacGrid, 2) - 1)
        TayMacMa = TayMacGrid(2 : end, 1, 1); % Free-stream Mach numbers used in the grid. dim(TayMacMa) = (size(TayMacGrid, 1) - 1) x 1
        TayMacSM = TayMacGrid(2 : end, 2 : end, :); % Shock angles and surface Mach numbers only (solutions to Taylor-Maccoll equations)
        [Xgrid, Ygrid] = meshgrid(TayMacOc, TayMacMa);
        % Use rational polynomial approximation for inverse P-M function (specific for Y = 1.4)
        vinf = (pi / 2) * (sqrt(6) - 1); % Maximum turning angle [rad]
        PrandltMeyerNu = @(M) sqrt((Y + 1) / (Y - 1)) * atan(sqrt((Y - 1) / (Y + 1) * (M^2 - 1))) - atan(sqrt(M^2 - 1)); % Prandyl-Meyer Angle [rad]
        PrandtlMeyerM2 = @(y) (1 + 1.3604 * y + 0.0962 * y^2 - 0.5127 * y^3) / (1 - 0.6722 * y - 0.3278 * y^2); % Gives Mach number for given ratio of PM angles y = (v / vinf)^(2/3) []
        % Define area ratio
        OmegaAreaRatio = @(M) ((2 / (Y + 1)) * (1 + (Y - 1) / 2 * M^2))^((1 / 2) * (Y + 1) / (Y - 1)) / M;
        numEl = 100;
        xCone = [];
        Pcone = [];
        % Find pressure on each cone
        for ii = 1 : numCones
            Oc = atan(mTangent(ii)); % Actual cone angle [rad]
            Os = interp2(Xgrid, Ygrid, TayMacSM(:, :, 1), Oc, Minf, 'linear'); % Shock angle at Minf
            Mc = interp2(Xgrid, Ygrid, TayMacSM(:, :, 2), Oc, Minf, 'linear'); % Mach number on cone
            % Find flow properties on the actual cone surface using isentropic & oblique shock equations
            Mninf = Minf * sin(Os); % Normal Mach number to shock []
            TtTinfRatio = 1 + (Y - 1) / 2 * Minf^2; % Isentropic flow, Tt = T0 through shock [K]
            Tt = TtTinfRatio * Tinf; % Total temperature (temperature on tip of nose cone) [K]
            Tcone = Tt / (1 + (Y - 1) / 2 * Mc^2); % Flow temperature along the cone surface [K]
            rhoSrhoinf = (Y + 1) * Mninf^2 / ((Y - 1) * Mninf^2 + 2); % Density ratio through shock []
            pSpinf = (2 * Y * Mninf^2 - (Y - 1)) / (Y + 1); % Pressure ratio through shock [Pa]
            ptinf = pinf * TtTinfRatio^(Y / (Y - 1)); % Total pressure in freestream [Pa]
            ptS = ptinf * rhoSrhoinf^(Y / (Y - 1)) * pSpinf^(1 / (1 - Y)); % Total pressure behind shock [Pa]
            pc(ii) = ptS * (1 + (Y - 1) / 2 * Mc^2)^(Y / (1 - Y)); % Equivalent tangent cone surface pressure [Pa]
            % Results for the current cone - Tt, Tcone, pc. Only first column is useful
            results(:, ii) = [pc(ii), Tt, Tcone]'; % Results for this cone. Tcone [K] and pc [Pa] along entire cone, Tt [K] on nose
            % Regions are split up between corners (1, 2, 3) separated by P-M expansion fans behind front cone
            xSection = linspace(xConnect(ii), xConnect(ii + 1), numEl); % Stations along x-axis in this section of cone [m]
            if (ii == 1)
                Psection = pc(1) * xSection;
            else
                delta2 = atan(mTangent(ii)) - atan(mTangent(ii + 1)); % delta - flow turning angle over previous corner [rad]
                d
            end
            xCone = [xCone, xSection]; % Stations along x-axis so far for entire nose cone [m]
            Pcone = [Pcone, Psection]; % Pressure distribution so far for entire nose cone [Pa]
        end
    end
end % Calculate aerodynamic coefficients in Ma, Re due to rocket geometry
function [Rlat, glat, g, Zg, T, p, rho, c, mu, nu, Cd, m, FT, FD, FG, FP, FN, v, q, Ma] = odetc(t, xx, args, seq, stage) 
    % Flight parameters - Stage 1 = S1 (Booster) & Stage 2 = S2 (Sustainer).
    Ar = args.rocketParams.Ar; % [S1 cross-sectional area, S2 cross-sectional area] [m2]
    mf = args.rocketParams.massFuelS; % [S1 fuel mass, S2 fuel mass] [kg]
    m0 = args.rocketParams.massInitS; % [Mass on ground, mass just after separation] [kg]
    hpc = args.rocketParams.hpc; % Hermite interpolating Polynomial Coefficients
    % Define auxiliary "functions" (constants at each step) for ODEs
    z = xx(3, :)'; % Altitude [m]
    v = sqrt(xx(4, :).^2 + xx(5, :).^2 + xx(6, :).^2)'; % Velocity magnitude [m/s]
    % Obtain grav. acceleration [m/s2], temperature [K], density [kg/m3], local speed of sound [m/s] at latitude and altitude
    [Rlat, glat, g, Zg, T, p, rho, c, mu, nu] = atmos(z + args.launchSiteParams.localElev, args.launchSiteParams.geodetLat, args.launchSiteParams.delT);
    q = (1 / 2) .* rho .* v.^2; % Dynamic pressure [Pa]
    Ma = v ./ c; % Mach number []
    % Recall csvs was defined {{MaCd{1}, tbpc{1}, tbFT{1}, tmdt{1}}, {MaCd{2}, tbpc{2}, tbFT{2}, tmdt{2}}}
    Cds = ppval(hpc{stage}{1}, Ma); % Un/powered drag coefficient []
%     Cd = coefficients(p, T, 1.00015, 0.01, args);
    numelt = numel(t); % Number of elements in (possible vector) t
    % Establish flight regime (staging/sequences) for appropriate forces FT, FD, FG
    switch seq
        case 1 % Stage 1 (on rail)
            [delr, delphi, Lcm, Ltower, L] = deal(args.initParams.delr, args.initParams.delphi, args.initParams.Lcm, args.launchSiteParams.towerSpan, args.rocketParams.allLength);
            phi = acos(z ./ sqrt(xx(1, :).^2' + xx(2, :).^2' + z.^2));
            rcmr = sqrt((xx(1, :)' + delr*cos(phi)).^2 + xx(2, :).^2' + (z - delr*sin(phi)).^2);
            Ltil = min(L(stage), max(0, Ltower - rcmr + Lcm));
            Smdtdt = min(0, ppval(hpc{stage}{5}, t));
            pe = nozzle(t, args, stage, p);
            % ~
            m = Smdtdt + m0(stage); % Mass [kg]
            Cd = Cds(2, :)'; % Powered drag coefficient []
            FT = ppval(hpc{stage}{3}, t) + (pe - p) * args.rocketParams.Ae(stage);
            FD = q .* Cd * Ar(1); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = max(0, FG * cosd(delphi) - FT); % Parallel rail [N]
            FN = FG .* (Ltil / L(stage)) * sind(delphi); % Normal rail [N]
        case 2 % Stage 1 (off rail)
            Smdtdt = min(0, ppval(hpc{stage}{5}, t));
            pe = nozzle(t, args, stage, p);
            % ~
            m = Smdtdt + m0(stage); % Mass [kg]
            Cd = Cds(2, :)'; % Powered drag coefficient []
            FT = max(0, ppval(hpc{stage}{3}, t) + (pe - p) * args.rocketParams.Ae(stage));
            FD = q .* Cd * Ar(stage); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 3 % Cruise with stage 1 attached (sequence 3)
            m = (m0(stage) - mf(stage)) * ones(numelt, 1); % S2 mass with all fuel and S1 attached [kg]
            Cd = Cds(1, :)'; % Unpowered drag coefficient []
            FT = zeros(numelt, 1); % Thrust [N]
            FD = q .* Cd * Ar(stage); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 4 % Cruise with stage 1 detached (sequence 3)
            m = m0(stage) * ones(numelt, 1); % S2 mass with all fuel and S1 attached [kg]
            Cd = Cds(1, :)'; % Unpowered drag coefficient []
            FT = zeros(numelt, 1); % Thrust [N]
            FD = q .* Cd * Ar(stage); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 5 % Stage 2 (sequence 5)
            [sepdly, igndly, tb] = deal(args.rocketParams.sepDelays, args.rocketParams.ignDelays, args.rocketParams.tb); % Burn time vector for stages 1, 2 [s]
            Smdtdt = min(0, ppval(hpc{stage}{5}, t - sepdly - igndly - tb(1)));
            pe = nozzle(t - sepdly - igndly - tb(1), args, stage, p);
            % ~
            m = Smdtdt + m0(stage); % Mass [kg]
            Cd = Cds(2, :)'; % Powered drag coefficient []
            FT = max(0, ppval(hpc{stage}{3}, t - sepdly - igndly - tb(1)) + (pe - p) * args.rocketParams.Ae(2));
            FD = q .* Cd * Ar(stage); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 6 % Cruising before apogee (sequence 6)
            m = (m0(stage) - mf(stage)) * ones(numelt, 1); % S2 with no fuel [kg]
            Cd = Cds(1, :)'; % Unpowered drag coefficient []
            FT = zeros(numelt, 1); % Thrust [N]
            FD = q .* Cd * Ar(stage); % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 7  % Descending with drogue chute (sequence 7)
            [AS2, Ap1, te6] = deal(args.rocketParams.Ar(stage), args.rocketParams.Ap(1), args.initParams.laste{1}(seq - 1));
            erf6 = erf((t - te6)/2);
            Cd = Cds(1, :)' + 0.97 * erf6; % Drogue parachute is deployed
            ADP = AS2 + Ap1 * erf6; % drogue descends at 115 ft/s
            % ~
            m = (m0(stage) - mf(stage)) * ones(numelt, 1); % S2 with no fuel [kg]
            FT = zeros(numelt, 1); % No thrust
            FD = q .* Cd .* ADP; % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        case 8 % Descending with parachute deployed (sequence 8)
            [AS2, Ap1, Ap2, te6, te7] = deal(args.rocketParams.Ar(2), args.rocketParams.Ap(1), args.rocketParams.Ap(2), args.initParams.laste{1}(seq - 2), args.initParams.laste{1}(seq - 1));
            [erf6, erf7] = deal(erf((t - te6)/2), erf((t - te7)/2));
            Cd = Cds(1, :)' + 0.97 * erf6 + 0.97 * erf7; % Main chute is deployed
            AP = AS2 + Ap1 * erf6 + Ap2 * erf7;
            % ~
            m = (m0(2) - mf(2)) * ones(numelt, 1); % S2 with no fuel [kg]
            FT = zeros(numelt, 1); % No thrust
            FD = q .* Cd .* AP; % Drag [N]
            FG = m .* g;  % Gravity [N]
            FP = zeros(numelt, 1); % Parallel rail [N]
            FN = zeros(numelt, 1); % Normal rail [N]
        otherwise
            disp("No sequence specified.")
    end
end  % Calculate forces
function dxxdt = odeval(t, xx, args, seq, stage) 
    % Reset variables for brevity
    StochXmax = args.stochasticParams.StochXmax;
    delphi = args.initParams.delphi;
    deltheta = args.initParams.deltheta;
    % Display time
    if (args.iWant.ODEsTimes), disp("t = " + t), end
    switch seq
        case 1 % No stochastics on the rail
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, m, FT, FD, FG, FP, FN, ~, ~, ~] = odetc(t, xx, args, seq, stage);
            [FXx, FXy, FXz] = deal(0); % No stochastics
            % Define direction sines and cosines in the form of geometric construction of launch rail
            vxv = sind(delphi) * cosd(deltheta);
            vyv = sind(delphi) * sind(deltheta);
            vzv = cosd(delphi);
        case {2, 3, 4, 5, 6, 7, 8}
            if (args.iWant.HaveStoch && args.iWant.RandStoch)
                StochXmax = 2 * args.stochasticParams.RandSXmax * (rand - 1);
            end
            [~, ~, ~, ~, ~, ~, rho, ~, ~, ~, ~, m, FT, FD, FG, FP, FN, v, ~, ~] = odetc(t, xx, args, seq, stage);
            % Define stochastic force FX
            [FXx, FXy, FXz] = deal(0); % Set stochastic forcing to 0 by default
            if (args.iWant.HaveStoch && abs(StochXmax) < 1 && ~(seq == 5 || seq == 6 || seq == 7))
                % Note: rand is a uniformly distributed random number on (0, 1)
                FXx = StochXmax * (2 * rand - 1) * (rho / args.launchSiteParams.rho0) * FG; % Random force in x
                FXy = StochXmax * (2 * rand - 1) * (rho / args.launchSiteParams.rho0) * FG; % Random force in y
                FXz = StochXmax * (2 * rand - 1) * (rho / args.launchSiteParams.rho0) * FG; % Random force in z
            elseif (args.iWant.HaveStoch && abs(StochXmax) >= 1)
                disp('Stochastic forcing should be perturbative. Adjust n so that |n| < 1.')
                return
            end
            % Apply a correction if necessary
            if (v == 0), v = v + eps; end % Perturb v slightly if delphi = 0 and odeinstep ~ 0 to avoid singularity
            % Define direction sines and cosines in the form of ratio of velocity components
            vxv = xx(4) / v;
            vyv = xx(5) / v;
            vzv = xx(6) / v;
        otherwise
            disp("No flight sequence specified")
    end
    % Define xx(1) = x, xx(2) = y, xx(3) = z, xx(4) = x', xx(5) = y', xx(6) = z'
    dxxdt = zeros(6, 1);
    dxxdt(1) = xx(4); % x' = x'
    dxxdt(2) = xx(5); % y' = y'
    dxxdt(3) = xx(6); % z' = z'
    dxxdt(4) = ((FT - FD) * vxv + FXx + (FP * sind(delphi) - FN * cosd(delphi)) * cosd(deltheta)) / m; % x'' = Fx / m
    dxxdt(5) = ((FT - FD) * vyv + FXy + (FP * sind(delphi) - FN * cosd(delphi)) * sind(deltheta)) / m; % x'' = Fy / m
    dxxdt(6) = ((FT - FD) * vzv - FG + FXz + FP * cosd(delphi) + FN * sind(delphi)) / m; % z'' = Fz / m
end  % Define xx' = f(t, xx) at time t
function [value, isterminal, direction] = odevents(t, xx, args, seq) 
    % See: https://www.mathworks.com/help/matlab/math/ode-event-location.html
    switch seq
        case 1
            % Create an event when the rail clears (angles no longer constant)
            [x, y, z] = deal(xx(1), xx(2), xx(3));
            r = sqrt(x^2 + y^2 + z^2);
            delr = args.initParams.delr;
            rcmr = sqrt((x + delr * z / r)^2 + y^2 + (z - delr * sqrt(1 - (z / r)^2))^2);
            value = rcmr - args.initParams.Lcm - args.launchSiteParams.towerSpan; % Rail is cleared
            isterminal = 1; % Stops the solver
            direction = 1; % Find increasing root of value
        case 2 % Create an event at the first (booster) burnout
            value = t - args.rocketParams.tb(1); % Burnout
            isterminal = 1; % Stops solver
            direction = 1; % Find increasing root
        case 3 % Create an event when the separation delay period ends if sepdly ~= 0
            value = t - (args.rocketParams.tb(1) + args.rocketParams.sepDelays); % Separation coast over
            isterminal = 1; % Stops solver
            direction = 1; % Find increasing root
        case 4 % Create an event when the ignition delay period ends if igndly ~= 0
            value = t - (args.rocketParams.tb(1) + args.rocketParams.sepDelays + args.rocketParams.ignDelays); % Ignition coast over
            isterminal = 1; % Stops solver
            direction = 1; % Find increasing root
        case 5 % Create an event at the second (sustainer) burnout
            value = t - (sum(args.rocketParams.tb) + args.rocketParams.sepDelays + args.rocketParams.ignDelays); % Burnout
            isterminal = 1; % Stops solver
            direction = 1; % Find increasing root
        case 6 % Create an event at apogee
            value = xx(6); % Apogee
            isterminal = 1; % Stops solver
            direction = -1; % Find decreasing root
        case 7 % Create an event when deploying the main chute
            value = xx(3) - args.rocketParams.deployAlt; % Main chute deployed
            isterminal = 1; % Stops solver
            direction = -1; % Find decreasing root
        case 8 % Create an event when touching down
            value = xx(3); % Impact
            isterminal = 1; % Stops solver - End of trajectory
            direction = -1; % Find decreasing root
    end
end  % Events that integrator tracks
function [dvx, dvy, dvz, dv] = delv(tau, t, acell, FXcell) 
    % tau is of the form [seq0, seq1] to specify which seqs over which to integrate
    seqs = min(tau) : max(tau);
    t = vertcat(t{seqs});
    amat = vertcat(acell{seqs});
    FXx = vertcat(FXcell{1}{seqs});
    FXy = vertcat(FXcell{2}{seqs});
    FXz = vertcat(FXcell{3}{seqs});
    FXvec = [FXx, FXy, FXz];
    if (any(any(FXvec)))
        nanmat = isnan(amat);
        [nanrows, ~] = find(nanmat);
        nanrows = unique(nanrows);
        amat(nanrows, :) = [];
        t(nanrows) = [];
    end
    dvvec = deal(trapz(t, amat));
    [dvx, dvy, dvz, dv] = deal(dvvec(1), dvvec(2), dvvec(3), norm(dvvec));
end % Estimate delta-V between sequences due to FT, FG, FD
function plots(Y, args) 
    % Begin plotting
    close all
    set(0,'DefaultFigureWindowStyle','docked')
    numelt = numel(vertcat(Y.t));
    blue = [0, 0.4470, 0.7410];
    orange = [0.8500, 0.3250, 0.0980];
    yellow = [0.9290, 0.6940, 0.1250]; 
    purple = [0.4940, 0.1840, 0.5560];
    green = [0.4660, 0.6740, 0.1880];
    lblue = [0.3010, 0.7450, 0.9330];
    red = [0.6350, 0.0780, 0.1840];
    mksz = 4;
    
    event = args.plottingtools.evseqs;
    SolnOnly = args.iWant.OmitCalcs;
    units = args.iWant.angleunit;
    
    t = vertcat(Y.t{:});
    x = vertcat(Y.x{:});
    y = vertcat(Y.y{:});
    z = vertcat(Y.z{:});
    vx = vertcat(Y.vx{:});
    vy = vertcat(Y.vy{:});
    vz = vertcat(Y.vz{:});
    
    delphi = args.initParams.delphi;
    deltheta = mod(args.initParams.deltheta, 360); % Ensure 0 <= deltheta < 360. Works with negative angles too. [deg]
    sepdly = args.rocketParams.sepDelays;
    igndly = args.rocketParams.ignDelays;
    % Find the direction of launch based on deltheta (0 = "East" and 90 = "North")
    if (deltheta == 11.25), direction = "E-ENEward";
    elseif (deltheta > 11.25 && deltheta < 33.75), direction = "East NEward";
    elseif (deltheta == 33.75), direction = "ENE-NEward";
    elseif (deltheta > 33.75 && deltheta < 56.25), direction = "NEward";
    elseif (deltheta == 56.25), direction = "NNE-NEward";
    elseif (deltheta > 56.25 && deltheta < 78.75), direction = "North NEward";
    elseif (deltheta == 78.75), direction = "N-NNEward";
    elseif (deltheta > 78.75 && deltheta < 101.25), direction = "Northward";
    elseif (deltheta == 101.25), direction = "NNW-Nward";
    elseif (deltheta > 101.25 && deltheta < 123.75), direction = "North NWward";
    elseif (deltheta == 123.75), direction = "NW-NNWward";
    elseif (deltheta > 123.75 && deltheta < 146.25), direction = "NWward";
    elseif (deltheta == 146.25), direction = "WNW-NWward";
    elseif (deltheta > 146.25 && deltheta < 168.75), direction = "West NWward";
    elseif (deltheta == 168.75), direction = "W-WNWward";
    elseif (deltheta > 168.75 && deltheta < 191.25), direction = "Westward";
    elseif (deltheta == 191.25), direction = "WSW-Wward";
    elseif (deltheta > 191.25 && deltheta < 213.75), direction = "West SWward";
    elseif (deltheta == 213.75), direction = "SW-WSWward";
    elseif (deltheta > 213.75 && deltheta < 236.25), direction = "SWward";
    elseif (deltheta == 236.25), direction = "SSW-SWward";
    elseif (deltheta > 236.25 && deltheta < 258.75), direction = "South SWward";
    elseif (deltheta == 258.75), direction = "S-SSWward";
    elseif (deltheta > 258.75 && deltheta < 281.25), direction = "Southward";
    elseif (deltheta == 281.25), direction = "SSE-Sward";
    elseif (deltheta > 281.25 && deltheta < 303.75), direction = "South SEward";
    elseif (deltheta == 303.75), direction = "SE-SSEward";
    elseif (deltheta > 303.75 && deltheta < 326.25), direction = "SEward";
    elseif (deltheta == 326.25), direction = "ESE-SEward";
    elseif (deltheta > 326.25 && deltheta < 348.75), direction = "East SEward";
    elseif (deltheta == 348.75), direction = "E-ESEward";
    else, direction = "Eastward";
    end
    
    if (SolnOnly) % Create plots with minimum information using only (t, xx)
        figure
        subplot(2, 1, 1) % Position (x, y, z)
        
        subplotStrT = sprintf('Cartesian State Representation to %s', event(numelt));
%         subplotStrB = sprintf('with $\\delta_\\phi = %2.2f^\\circ, \\ \\delta_\\theta = %2.2f^\\circ,$ and $%2.2f$ [s] Delay', delphi, deltheta, delay);
        subplotStrB = sprintf('with %s $%2.2f^\\circ$ Launch Angle and $%2.2f$ [s] Sep. and $%2.2f$ [s] Ign. Delays', direction, delphi, sepdly, igndly);
        sgtitle({subplotStrT, subplotStrB}, 'interpreter', 'latex')
        title("Position")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, z / 1000, '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, x / 1000, '.', 'markersize', mksz, 'color', orange)
        plot(t, y / 1000, '.', 'markersize', mksz, 'color', red)
        axes()
        yyaxis left
        ylabel("$z$ [$km$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$x, \\ y$ [$km$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legend(legl, {'$x$ [$km$]', '$y$ [$km$]', '$z$ [$km$]'}, 'interpreter', 'latex');
        % v^v^v^v
        subplot(2, 1, 2) % Velocity (vx, vy, vz)
        title("Velocity")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vz, '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, vx, '.', 'markersize', mksz, 'color', orange)
        plot(t, vy, '.', 'markersize', mksz, 'color', red)
        axes()
        yyaxis left
        ylabel("$v_z$ [$m/s$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$v_x, \\ v_y$ [$m/s$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legend(legl, {'$v_x$ [$m/s$]', '$v_y$ [$m/s$]', '$v_z$ [$m/s$]'}, 'interpreter', 'latex');
        
        % Plot z vs r*sin(phi)
        figure
        title('Geometric Trajectory')
        xlabel('Horizontal Displacement $r \sin\phi$ [$km$]', 'interpreter', 'latex')
        ylabel('Altitude $z$ [$km$]', 'interpreter', 'latex')
        grid on; hold on
        rxy = sqrt(x.^2 + y.^2); % Perform minor calculation to find horz. distance
        plot(rxy / 1000, z / 1000, '.', 'markersize', mksz, 'color', blue)
        axis equal
        
        % Plot v vs r
        figure
        title('Radial Displacement and Speed')
        xlabel("Speed $v$ [$m/s$]", 'interpreter', 'latex')
        ylabel("Radial Displacement $r$ [$km$]", 'interpreter', 'latex')
        grid on; hold on
        r = sqrt(x.^2 + y.^2 + z.^2); % Perform minor calculation to find distance
        v = sqrt(vx.^2 + vy.^2 + vz.^2); % Perform minor calculation to find speed
        plot(v, r/1000, '.', 'markersize', mksz, 'color', blue)
        
        
    else % Use data obtained from the solutions (with calculations)
        r = vertcat(Y.r{:});
        v = vertcat(Y.v{:});
        q = vertcat(Y.q{:});
        FD = vertcat(Y.FD{:});
    % Plot the full state (pos, vel, accel) together as one figure
        figure
        subplot(3, 1, 1) % Position (x, y, z)
        subplotStrT = sprintf('Cartesian State Representation to %s', event(numelt));
%         subplotStrB = sprintf('with $\\delta_\\phi = %2.2f^\\circ, \\ \\delta_\\theta = %2.2f^\\circ,$ and $%2.2f$ [s] Delay', delphi, deltheta, delay);
        subplotStrB = sprintf('with %s $%2.2f^\\circ$ Launch Angle and $%2.2f$ [s] Sep., $%2.2f$ [s] Ign. Delays', direction, delphi, sepdly, igndly);
        sgtitle({subplotStrT, subplotStrB}, 'interpreter', 'latex')
        title("Position")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, z / 1000, '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, x / 1000, '.', 'markersize', mksz, 'color', orange)
        plot(t, y / 1000, '.', 'markersize', mksz, 'color', red)
        axes()
        yyaxis left
        ylabel("$z$ [$km$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$x, \\ y$ [$km$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legend(legl, {'$x$ [$km$]', '$y$ [$km$]', '$z$ [$km$]'}, 'interpreter', 'latex');
        % v^v^v^v
        subplot(3, 1, 2) % Velocity (vx, vy, vz)
        title("Velocity")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vz, '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, vx, '.', 'markersize', mksz, 'color', orange)
        plot(t, vy, '.', 'markersize', mksz, 'color', red)
        axes()
        yyaxis left
        ylabel("$v_z$ [$m/s$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$v_x, \\ v_y$ [$m/s$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legend(legl, {'$v_x$ [$m/s$]', '$v_y$ [$m/s$]', '$v_z$ [$m/s$]'}, 'interpreter', 'latex');
        % v^v^v^v
        subplot(3, 1, 3) % Acceleration (ax, ay, az)
        title("Acceleration")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vertcat(Y.az{:}), '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, vertcat(Y.ax{:}), '.', 'markersize', mksz, 'color', orange)
        plot(t, vertcat(Y.ay{:}), '.', 'markersize', mksz, 'color', red)
        axes()
        yyaxis left
        ylabel("$a_z$ [$m/s^2$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$a_x, \\ a_y$ [$m/s^2$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legend(legl, {'$a_x$ [$m/s^2$]', '$a_y$ [$m/s^2$]', '$a_z$ [$m/s^2$]'}, 'interpreter', 'latex');
        clear legl


        % Plot the forces all together - No need to state the seq, delphi, deltheta anymore
        if (args.iWant.HaveStoch)
            figure
            hold on; grid on
            title("Forces")
            xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
            yyaxis left
            plot(t, vertcat(Y.FT{:}), '.', 'markersize', mksz, 'color', blue)
            plot(t, FD, '.', 'markersize', mksz, 'color', lblue)
            plot(t, vertcat(Y.FG{:}), '.', 'markersize', mksz, 'color', green)
            yyaxis right
            plot(t, vertcat(Y.FXx{:}), '.', 'markersize', mksz, 'color', orange)
            plot(t, vertcat(Y.FXy{:}), '.', 'markersize', mksz, 'color', yellow)
            plot(t, vertcat(Y.FXz{:}), '.', 'markersize', mksz, 'color', red)
            axes()
            yyaxis left
            ylabel("Body Forces $F_T,\ F_D,\ F_G$ [$N$]", 'interpreter', 'latex')
            set(gca, 'YTickLabel', get(gca, 'YTick'));
            yyaxis right
            ylabel("Stochastic Forces $F_{\chi_x}),\ F_{\chi_y}),\ F_{\chi_z})$ [$N$]", 'interpreter', 'latex')
            legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
            legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', lblue);
            legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', green);
            legl(4) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
            legl(5) = plot(nan, nan, '-', 'markersize', mksz, 'color', yellow);
            legl(6) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
            legend(legl, {'$F_T$', '$F_D$', '$F_G$', '$F_{\chi_x})$', '$F_{\chi_y})$', '$F_{\chi_z})$'}, 'interpreter', 'latex','NumColumns', 2)
            clear legl
        else
            figure
            hold on; grid on
            title("Forces")
            xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
            for ii = 1:numelt
                yyaxis left
                plot(t, vertcat(Y.FT{:}), '.', 'markersize', mksz, 'color', blue)
                yyaxis right
                plot(t, FD, '.', 'markersize', mksz, 'color', orange)
                plot(t, vertcat(Y.FG{:}), '.', 'markersize', mksz, 'color', red)
            end
            yyaxis left; yliml = get(gca,'Ylim');
            set(gca,'Ylim',[0, yliml(2)])
            yyaxis right; ylimr = get(gca,'Ylim');
            set(gca,'Ylim',[0, ylimr(2)])
            clear yliml ylimr
            yyaxis left
            ylabel("Accelerating Forces $F_T$ [$N$]", 'interpreter', 'latex')
            set(gca, 'YTickLabel', get(gca, 'YTick'));
            yyaxis right
            ylabel("Decelerating Forces $F_D,\ F_G$ [$N$]", 'interpreter', 'latex')
            legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
            legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
            legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', red);
            legend(legl, {'$F_T$ [$N$]', '$F_D$ [$N$]', '$F_G$ [$N$]'}, 'interpreter', 'latex')
            clear legl
        end



        % Plot the angles all together
        figure
        hold on; grid on
        title("Angles in Origin and Body Frames")
        xlabel('Time $t$ [$s$]', 'interpreter', 'latex')
        padd = 28.648/2; % Padding in degrees
        padr = 0.5/2; % Padding in radians
        if (units == "d" || units == "deg" || units == "degree" || units == "degrees")
            plot(t, vertcat(Y.thetad{:}), '.', 'markersize', mksz, 'color', blue)
            plot(t, vertcat(Y.phid{:}), '.', 'markersize', mksz, 'color', orange)
            plot(t, vertcat(Y.thetapd{:}), '.', 'markersize', mksz, 'color', yellow)
            plot(t, vertcat(Y.phipd{:}), '.', 'markersize', mksz, 'color', purple)
            plot(t, vertcat(Y.psid{:}), '.', 'markersize', mksz, 'color', green)
            ylabel('Angles $\theta, \phi, \theta'', \phi'', \psi$ [deg]', 'interpreter', 'latex')
            ylim([-180 - padd, 180 + padd]);
            set(gca,'YTick',-180:45:180) 
            set(gca,'YTickLabel', {'-180', '-135', '-90', '-45', '0', '45', '90', '135', '180'})
        elseif (units == "r" || units == "rad" || units == "radian" || units == "radians")
            plot(t, vertcat(Y.theta{:}), '.', 'markersize', mksz, 'color', blue)
            plot(t, vertcat(Y.phi{:}), '.', 'markersize', mksz, 'color', orange)
            plot(t, vertcat(Y.thetap{:}), '.', 'markersize', mksz, 'color', yellow)
            plot(t, vertcat(Y.phip{:}), '.', 'markersize', mksz, 'color', purple)
            plot(t, vertcat(Y.psi{:}), '.', 'markersize', mksz, 'color', green)
            ylabel('Angles $\theta, \phi, \theta'', \phi'', \psi$ [rad]', 'interpreter', 'latex')
            ylim([-padr, pi + padr]);
            set(gca,'YTick',0:pi/8:pi) 
            set(gca,'YTickLabel', {'0', '{\fontsize{11}\pi}/8', ...
            '{\fontsize{11}\pi}/4', '{\fontsize{9}3}{\fontsize{11}\pi}/8', '{\fontsize{11}\pi}/2', ...
            '{\fontsize{9}5}{\fontsize{11}\pi}/8', '{\fontsize{9}3}{\fontsize{11}\pi}/4', ...
            '{\fontsize{9}7}{\fontsize{11}\pi}/8', '{\fontsize{11}\pi}'})
        else
            disp("Check units in Switches is any of d, deg, degree, degrees or r, rad, radian, radians.")
        end
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', yellow);
        legl(4) = plot(nan, nan, '-', 'markersize', mksz, 'color', purple);
        legl(5) = plot(nan, nan, '-', 'markersize', mksz, 'color', green);
        legend(legl, {'$\theta$', '$\phi$', '$\theta''$', '$\phi''$', '$\psi$'}, 'interpreter', 'latex')
        clear legl



        % Plot dynamic angles in flight
        figure
        grid on; hold on
        title("Thrust and Drag Directionality Functions")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        plot(t, vertcat(Y.vxv{:}), '.', 'markersize', mksz, 'color', blue)
        plot(t, vertcat(Y.vyv{:}), '.', 'markersize', mksz, 'color', orange)
        plot(t, vertcat(Y.vzv{:}), '.', 'markersize', mksz, 'color', yellow)
        ylabel("Modulating Velocity Components $v_{x,y,z}/v$", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legl(3) = plot(nan, nan, '-', 'markersize', mksz, 'color', yellow);
    %     legend(legl, {'$v_x/v$', '$v_y/v$', '$v_z/v$'}), 'interpreter', 'latex')
        legend(legl, ["$\sin \phi' \cos \theta'$", "$\sin \phi' \sin \theta'$", "$\cos \phi'$"], 'interpreter', 'latex')
        clear legl

        % Plot flight parameters together
        figure
        subplot(3, 1, 1) % Mass
        sgtitle('Flight Parameters')
        title("Mass")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, vertcat(Y.m{:}), '.', 'markersize', mksz, 'color', blue)
        ylimm = ylim;
        ylim([0, ylimm(2)]);
        clear ylimm
        ylabel("$m$ [$kg$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(3, 1, 2) % Re and Cd
        title("Reynolds Number and Drag Coefficient")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vertcat(Y.Re{:}), '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, vertcat(Y.Cd{:}), '.', 'markersize', mksz, 'color', orange)
        yyaxis left
        ylabel("$Re$ [$\,$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$C_d$ [$\,$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legend(legl, {'$Re$ [$\,$]', '$C_d$ [$\,$]'}, 'interpreter', 'latex');
        % v^v^v^v
        subplot(3, 1, 3) % Ma and q
        title("Mach Number and Dynamic Pressure")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vertcat(Y.Ma{:}), '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, q / 1000, '.', 'markersize', mksz, 'color', orange)
        axes()
        yyaxis left
        ylabel("$Ma$ [$\,$]", 'interpreter', 'latex')
        yyaxis right
        ylabel("$q$ [$kPa$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legend(legl, {'$Ma$ [$\,$]', '$q$ [$kPa$]'}, 'interpreter', 'latex');
        clear legl



        % Plot the atmosphere
        figure
        subplot(3, 2, 1)
        sgtitle('Local Atmospheric Properties of Air During Flight')
        title("Temperature")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, vertcat(Y.T{:}), '.', 'markersize', mksz, 'color', blue)
        ylabel("$T$ [$K$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(3, 2, 2)
        title("Speed of Sound")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, vertcat(Y.c{:}), '.', 'markersize', mksz, 'color', blue)
        ylabel("$c$ [$m/s$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(3, 2, 3)
        title("Pressure")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, vertcat(Y.p{:})/1000, '.', 'markersize', mksz, 'color', blue)
        ylabel("$p$ [$kPa$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(3, 2, 4)
        title("Viscosity")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        yyaxis left
        plot(t, vertcat(Y.mu{:}), '.', 'markersize', mksz, 'color', blue)
        yyaxis right
        plot(t, vertcat(Y.nu{:}), '.', 'markersize', mksz, 'color', orange)
        yyaxis left
        ylabel("Dynamic $\mu$ [kg/m $\cdot$ s]", 'interpreter', 'latex')
        yyaxis right
        ylabel("Kinematic $\nu$ [$m^2/s$]", 'interpreter', 'latex')
        legl(1) = plot(nan, nan, '-', 'markersize', mksz, 'color', blue);
        legl(2) = plot(nan, nan, '-', 'markersize', mksz, 'color', orange);
        legend(legl, {'Dynamic $\mu$ [$kg/m\cdot s$]', 'Kinematic $\nu$ [$m^2/s$]'}, 'interpreter', 'latex')
        % v^v^v^v
        subplot(3, 2, [5, 6])
        title("Density")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, vertcat(Y.rho{:}), '.', 'markersize', mksz, 'color', blue)
        ylabel("$\rho$ [$kg/m^3$]", 'interpreter', 'latex')



        % More customized plots

        % Plot time and velocity against force (left) and coefficient (right) of drag
        figure
        subplot(2, 1, 1)
        sgtitle('Drag')
        title("in Time")
        xlabel("Time $t$ [$s$]", 'interpreter', 'latex')
        ylabel("$F_D$ [$kN$]", 'interpreter', 'latex')
        hold on; grid on
        plot(t, FD/1000, '.', 'markersize', mksz, 'color', blue)
        % v^v^v^v
        subplot(2, 1, 2)
        title("in Speed")
        xlabel("Speed $v$ [$m/s$]", 'interpreter', 'latex')
        ylabel("$F_D$ [$kN$]", 'interpreter', 'latex')
        hold on; grid on
        plot(v, FD/1000, '.', 'markersize', mksz, 'color', blue)

        % Plot dynamic pressure against time, velocity, and altitude
        figure
        subplot(1, 3, 1) % Mass
        sgtitle('Dynamic Pressure')
        title("in Time")
        ylabel("Time $t$ [$s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(q/1000, t, '.', 'markersize', mksz, 'color', blue)
        xlabel("$q$ [$kPa$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(1, 3, 2) % Re and Cd
        title("in Speed")
        ylabel("Speed $v$ [$m/s$]", 'interpreter', 'latex')
        hold on; grid on
        plot(q/1000, v, '.', 'markersize', mksz, 'color', blue)
        xlabel("$q$ [$kPa$]", 'interpreter', 'latex')
        % v^v^v^v
        subplot(1, 3, 3) % Ma and q
        title("in Altitude")
        ylabel("Altitude $z$ [$km$]", 'interpreter', 'latex')
        hold on; grid on
        plot(q/1000, z/1000, '.', 'markersize', mksz, 'color', blue)
        xlabel("$q$ [$kPa$]", 'interpreter', 'latex')


        % Plot z vs r*sin(phi)
        figure
        title('Geometric Trajectory')
        xlabel('Horizontal Displacement $r \sin\phi$ [$km$]', 'interpreter', 'latex')
        ylabel('Altitude $z$ [$km$]', 'interpreter', 'latex')
        grid on; hold on
        plot(vertcat(Y.rxy{:}) / 1000, z/1000, '.', 'markersize', mksz, 'color', blue)
        axis equal


        % Plot v vs r
        figure
        title('Radial Displacement and Speed')
        xlabel("Speed $v$ [$m/s$]", 'interpreter', 'latex')
        ylabel("Radial Displacement $r$ [$km$]", 'interpreter', 'latex')
        grid on; hold on
        plot(v, r/1000, '.', 'markersize', mksz, 'color', blue)
        
        
        % Plot v vs z
        figure
        title('Altitude and Speed')
        xlabel("Speed $v$ [$m/s$]", 'interpreter', 'latex')
        ylabel("Altitude $z$ [$km$]", 'interpreter', 'latex')
        grid on; hold on
        plot(v, z/1000, '.', 'markersize', mksz, 'color', blue)
    end
        

    % Start window at first figure
    figure(1)
    % Reset default window docking property
    set(0,'DefaultFigureWindowStyle','normal')
        
    
    function axes()
        yyaxis right; ylimr = get(gca,'Ylim'); ratior = ylimr(1)/ylimr(2);
        yyaxis left; yliml = get(gca,'Ylim'); ratiol = yliml(1)/yliml(2);
        if (ratior ~= 0)
            if (yliml(2) * ratior < yliml(1))
                set(gca,'Ylim',[yliml(2)*ratior, yliml(2)])
            else
                set(gca,'Ylim',[yliml(1), yliml(1)/ratior])
            end
        elseif (ratiol ~= 0)
            yyaxis right
            if (ylimr(2) * ratiol < ylimr(1))
                set(gca,'Ylim',[ylimr(2)*ratiol, ylimr(2)])
            else
                set(gca,'Ylim',[ylimr(1), ylimr(1)/ratiol])
            end
        end
    end  % Fix plot axes to align at 0
end  % Plots function
