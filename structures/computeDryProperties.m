function rocket = computeDryMassProperties(rocket)
% 
% Matt Werner (m.werner@vt.edu) - Jan 14, 2021
% 
% Determine the dry mass of the rocket by calculating the (constant-
% density) material's volume, multiplying through by the corresponding
% material densities, and finally summing masses over the components. 
% 
% Also calculate the components' center of mass positions relative to 
% the coordinate system whose origin is placed along each component's
% centerline at the forwardmost location. The x-axis of this coordinate
% system points towards the body's center of mass and the z-axis points
% towards where the launch rail would be if the rocket were set for launch.
% The overall center of mass is obtained via the simple weighted sum
% definition, but having all positions taken with respect to the tip of the
% nosecone. Thus, components of the rocket are assumed to be ordered in the
% following fashion: nosecone, frustum, body, frustum, body, frustum, ...
% where the pattern repeats for as many stages as indicated.
% 
% Also calculate the moment of inertia matrix relative to this same
% coordinate frame used for the center of mass by finding the individiual
% inertia matrices of the components and repeated use of the 3D-parallel
% axis theorem.
% 
% That is, calculate:
% 1. The total dry mass (mass of the vehicle without fuel/motor) for
%    each stage
% 2. The position of the center of (dry) mass relative to the nosecone's
%    tip for each stage. The axis tracking the position of the center of
%    mass is the body-fixed principal axis frame shifted from the center of
%    mass (where the principal frame is defined) to the nosecone's tip.
%    This principal frame is oriented such that its z-axis points upwards
%    along the body's line of symmetry out through the nosecone's tip and
%    its x-axis (which is normal to the z-axis) points towards the rail
%    during the time prior to launch.
% 3. The moment of inertia matrix of the dry mass vehicle for each stage
%    with respect to this principal-frame axis positioned on the nosecone's
%    tip as described previously. By choosing the principal frame as the
%    body-fixed axis, the products of inertia for the flight vehicle are
%    zero.
%    Note: This moment of inertia matrix is NOT(!) taken with respect to
%    the vehicle's center of mass.
% 
% NOTE: There are 2 different orientations for coordinate systems here. The
%       principal axis frame is defined with its z-axis pointing up through
%       the nose cone and x-axis pointing towards where the launch rail
%       would be if the rocket were staged for launch, 
%       BUT
%       the coordinate frame used to find individual component centers of
%       mass and inertia matrices has its x-axis pointed into the rocket
%       (opposite direction of the principal frame z-axis) and its z-axis
%       pointing towards where the rail would be. Thus, a (simple)
%       coordinate transformation is required for both the official center
%       of mass and inertia matrix for them to be reported with respect to
%       the principal axis frame.
% 
%    Inputs:
% 
%            rocket - The input details of the rocket in standard SI units
%                     concerning the nose cone (.nosecone), body (.body),
%                     etc. These parameters are used in determining the
%                     rocket's structural properties (mass, center of mass,
%                     and moment of inertia matrix).
%                     Size: ? (structure)
%                     Units: SI
% 
%    Outputs:
% 
%            rocket - The newly updated structure containing all previous
%                     fields of the provided rocket plus the newly-found
%                     total mass, center of mass, and inertia matrix.
%                     Size: ? (structure)
%                     Units: SI
% 

% Required transformation to flip the BOR coordinate frame into the
% body-fixed principal frame
Tbody_bor = getTransformationBOR2Body;

%% Nosecone
% There is only a single nosecone (independent of the number of stages
% within the rocket)
[rocket.nosecone.profile, ...
    rocket.nosecone.mass, ...
    rocket.nosecone.CoM_bor, ...
    rocket.nosecone.atCoM.IMoI_bor] ...
    = computeNoseconeProperties(rocket.nosecone);

%% Frustum & body components
% Cycle through each frustum and body piece starting from the last stage
% (the pieces connected to the nosecone) to the first stage (the pieces
% connected to the booster)
for stage = rocket.stages:-1:1
    % Find the structural properties of this stage's body
    [rocket.body.profile{stage, 1}, ...
        rocket.body.mass(stage, 1), ...
        rocket.body.CoM_bor{stage, 1}, ...
        rocket.body.atCoM.IMoI_bor{stage, 1}] ...
        = computeCylinderProperties(rocket.body, stage);
    
    % Find the structural properties of this stage's frustum
    [rocket.frustum.profile{stage, 1}, ...
        rocket.frustum.mass(stage, 1), ...
        rocket.frustum.CoM_bor{stage, 1}, ...
        rocket.frustum.atCoM.IMoI_bor{stage, 1}] ...
        = computeFrustumProperties(rocket.frustum, stage);
end


%% Total properties
% Find the location of every component's mass center with respect to the
% nosecone's BOR coordinate system defined at the nose cone's tip with its
% x-axis pointing pointing along the nosecone's longitudinal axis towards
% its center of mass. Note that the entire dry mass rocket is aligned on
% the nosecone's centerline, so determining these positions occurs only on
% the x-axis. Also calculate each component's mass moment of inertia matrix
% with respect to the nosecone's BOR coordinate frame.
pCurrent_bor = [rocket.nosecone.length; 0; 0];
for stage = rocket.stages:-1:1
    % Every stage is paired with a frustum and body, so both may be
    % evaluated within this stage
    %
    % Find how far away this frustum's center of mass is from the nosecone
    % tip
    CoM_bor = pCurrent_bor + rocket.frustum.CoM_bor{stage, 1};
    rocket.frustum.atTip.CoM_bor{stage, 1} = CoM_bor;
    % Translate this frustum's mass moment of inertia matrix to the
    % nosecone tip
    rocket.frustum.atTip.IMoI_bor{stage, 1} = ...
        rocket.frustum.atCoM.IMoI_bor{stage, 1} + ...
        rocket.frustum.mass(stage, 1) * (CoM_bor'*CoM_bor*eye(3) - ...
                                         CoM_bor*CoM_bor');
    % Update the current position from the nosecone's tip by the frustum's
    % length to get to the body's coordinate origin
    pCurrent_bor = pCurrent_bor + [rocket.frustum.length(stage, 1); 0; 0];
    % Find how far away this cylinder/body's center of mass is from the
    % nosecone tip
    rocket.body.atTip.CoM_bor{stage, 1} = ...
        pCurrent_bor + rocket.body.CoM_bor{stage, 1};
    % Update the current position from from the nosecone's tip by the
    % cylinder/body's length to get to the next frustum's coordinate origin
    pCurrent_bor = pCurrent_bor + [rocket.body.length(stage, 1); 0; 0];
end

% Prepare for iterative procedure for finding the total mass, center of
% mass, and mass moment of inertia of the dry rocket
nextStageMass = rocket.nosecone.mass;
nextStageCoM_bor = rocket.nosecone.CoM_bor;
nextStagesumMassCoM_bor = nextStageMass*nextStageCoM_bor;

% Add these preparations as the final stage in case the nosecone is
% requested to separate from the final firing stage
rocket.dry.mass(rocket.stages+1, 1) = nextStageMass;
rocket.dry.CoM_bor{rocket.stages+1, 1} = nextStageCoM_bor;

% Assemble all pieces together and find the dry-mass, dry center of mass,
% and the dry mass moment of inertia for each stage of the vehicle
for stage = rocket.stages:-1:1
    %% Total mass
    % Addition of mass due to this interstage and casing
    massOfThisInterstage = rocket.frustum.mass(stage, 1);
    massOfThisCasing = rocket.body.mass(stage, 1);
    % Total mass of this stage
    massOfThisStage = nextStageMass + ...
        massOfThisInterstage + ...
        massOfThisCasing;
    % Mark the total mass of this stage as an official calculated value
    % for the entire rocket
    rocket.dry.mass(stage, 1) = massOfThisStage;
    
    %% Total center of mass
    % Weighted sum of masses and centers of mass for this interstage and
    % casing relative to the nosecone's BOR coordinate frame
    CoMofThisInterstage = rocket.frustum.atTip.CoM_bor{stage, 1};
    CoMofThisCasing = rocket.body.atTip.CoM_bor{stage, 1};
    sumMCoMofThisStage = nextStagesumMassCoM_bor + ...
        massOfThisInterstage*CoMofThisInterstage + ...
        massOfThisCasing*CoMofThisCasing;
    % Compute the position of the center of mass for this stage relative to
    % the nosecone BOR coordinate frame
    CoMofThisStage = sumMCoMofThisStage/massOfThisStage;
    % Mark the position of the center of mass as an official calculated
    % value for the entire rocket
    rocket.dry.CoM_bor{stage, 1} = CoMofThisStage;
    
    %% Total mass moment of inertia (@nosecone tip, BOR frame)
    
    
    %% Update next values
    % Update the next stage mass
    nextStageMass = massOfThisStage;
    % Update the next stage center of mass
    nextStagesumMassCoM_bor = sumMCoMofThisStage;
end


%% Compute the total dry mass of the vehicle