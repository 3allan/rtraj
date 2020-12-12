% One-time run script that calculates some additional parameters from the
% inputs of rtraj

%% Create flags
flag_ProvidedMassFlowRate = 1; % User input includes nonlinear dmdt = f(t)
flag_ProvidedChamberPressure = 1; % User input includes nonlinear pc = f(t)

%% Propulsion
% Correct thrust profile data from .csv file by removing excess data for
% times indicated after the burn times according to "burnTimes" and then
% ensure that the units are SI. Also ensure that the mass flow rate is
% negative.
for stage = stageNums'
    % Adjust thrust profile (necessary)
    FTProfile{stage}(FTProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
    FTProfile{stage} = convertxUnits(FTProfile{stage}, FTUnits__{stage}, ["s", "N"]);
    
    % Adjust mass flow rate profile
    if (~all(isnan(MFProfile{stage})))
        MFProfile{stage}(MFProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
        MFProfile{stage} = convertxUnits(MFProfile{stage}, MFUnits__{stage}, ["s", "kg/s"]);
        % Ensure that the mass flow rate is negative
        MFProfile{stage}(:, 2) = -abs(MFProfile{stage}(:, 2));
    else
        % If not provided, then use linear variation
        MFProfile{stage} = -massMotor(stage)/burnTimes(stage);
        flag_ProvidedMassFlowRate = 0;
    end
    
    % Adjust chamber pressure profile
    if (~all(isnan(PCProfile{stage})))
        PCProfile{stage}(PCProfile{stage}(:, 1) > burnTimes(stage, 1), :) = [];
        PCProfile{stage} = convertxUnits(PCProfile{stage}, PCUnits__{stage}, ["s", "Pa"]);
    else
        flag_ProvidedChamberPressure = 0;
    end
end

% Convert given diameters (easily measurable, but not used directly in the
% process of obtaining the solution) to circular areas (not so easily
% measurable, but used directly in the process of obtaining the solution)
areaRefer = pi*diaOuter_.^2/4; % Reference area for drag [m2]
areaThrt_ = pi*diaThroat.^2/4; % Area of the throat for propulsion [m2]
areaExit_ = pi*diaExit__.^2/4; % Area of nozzle exit for adjusting thrust [m2]
areaPara_ = pi*diaFlatDM.^2/8; % Inflated parachute area for drag on descent [m2]

% Obtain (H)ermite (P)olynomial (C)oefficients (HPC) for fast interpolation of
% the thrust curve, mass flow rate, and chamber pressure (if provided) per
% stage ('pchip')
% Key:
% stage, 1 - Thrust profile (Always required, so it goes first)
% stage, 2 - Mass flow rate (Optional, but still always there)
% stage, 3 - Chamber pressure (NaN if no chamber pressure is indicated)
HPC = cell(max(stageNums), 3);
for stage = stageNums'
    % Obtain (P)iecewise (C)ubic (H)ermite (I)nterpolating (P)olynomial
    % (PCHIP) coefficients 
    %
    % PCHIP the thrust profile
    HPC{stage, 1} = pchip(FTProfile{stage}(:, 1), FTProfile{stage}(:, 2));
    %
    % PCHIP the mass flow rate
    HPC{stage, 2} = pchip(MFProfile{stage}(:, 1), MFProfile{stage}(:, 2));
    %
    % PCHIP the chamber pressure, if provided
    if (flag_ProvidedChamberPressure)
        HPC{stage, 3} = pchip(PCProfile{stage}(:, 1), PCProfile{stage}(:, 2));
    else
        % Otherwise, give it NaN so that attempting to interpolate upon
        % HPC{..., 3} will result in a Matlab error with ppval
        HPC{stage, 3} = NaN;
    end
end
% Rename for inserting into the table 'pars'
FThpc = HPC(:, 1);
MFhpc = HPC(:, 2);
PChpc = HPC(:, 3);

%% Earth
% Convert indicates launch longitude and latitude to radians
LonLaunchrad = deg2rad(LonLaunch);
LatLaunchrad = deg2rad(LatLaunch);


% Compute the Earth's rotation rate vector with respect to the ENV frame's
% orientation centered at the Earth's core
w_ecef = [0; 0; w];
Tenv_ecf = getTransformationECF2ENVCoordinates(LonLaunchrad, LatLaunchrad);
w_env = Tenv_ecf*w_ecef;
% Rename by transposing the above results to pass along with pars
[w_ecef_rowvec, w_env_rowvec] = deal(w_ecef', w_env');
% Obtain the inverse rotation for future's sake
Tecf_env = Tenv_ecf';


% Determine the ECEF position of the launch site based on the given inputs
RNLaunch = computePrimeVerticalRadius(Req, e, LatLaunchrad, "geodetic");
RNLaunch_plus_hLaunch = RNLaunch + ElevationGPS;
xLaunch_ecf = RNLaunch_plus_hLaunch*cos(LatLaunchrad)*cos(LonLaunchrad);
yLaunch_ecf = RNLaunch_plus_hLaunch*cos(LatLaunchrad)*sin(LonLaunchrad);
zLaunch_ecf = ((1 - e^2)*RNLaunch + ElevationGPS)*sin(LatLaunchrad);
rLaunch_ecef = [xLaunch_ecf; yLaunch_ecf; zLaunch_ecf];
rLaunch_ecef_rowvec = rLaunch_ecef';
clear RNplush xLaunch_ecf yLaunch_ecf zLaunch_ecf