function rocket = computeDryProperties(rocket)
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
% following fashion: nosecone, shoulder, cylinder, shoulder, cylinder, shoulder, ...
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
%                     concerning the nose cone (.nosecone), cylinder (.cylinder),
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
    rocket.nosecone.IMoIatCoM_bor] ...
    = computeNoseconeProperties(rocket.nosecone);

%% Shoulder & cylinder components
% Cycle through each shoulder and cylinder piece starting from the last stage
% (the pieces connected to the nosecone) to the first stage (the pieces
% connected to the booster)
tipToOrigin = [rocket.nosecone.length; 0; 0];
for stage = rocket.stages:-1:1
    % ----------------------- SHOULDER -----------------------
    % Find the structural properties of this stage's shoulder
    [rocket.shoulder.profile{stage, 1}, ...
        rocket.shoulder.mass(stage, 1), ...
        rocket.shoulder.CoM_bor{stage, 1}, ...
        rocket.shoulder.IMoIatCoM_bor{stage, 1}] ...
        = computeFrustumProperties(rocket.shoulder, stage);
    
    % Find how far away the shoulder's center of mass is from the nosecone's
    % tip & update the distance from the nosecone's tip to the shoulder's end
    [tipToCoM, tipToOrigin] = ...
        distanceFromNoseconeTip(tipToOrigin, ...
        rocket.shoulder.CoM_bor{stage, 1}, rocket.shoulder.length(stage, 1));
    rocket.shoulder.tipToCoM_bor{stage, 1} = tipToCoM;


    % ----------------------- CYLINDER -----------------------
    % Find the structural properties of this stage's casing
    [rocket.cylinder.profile{stage, 1}, ...
        rocket.cylinder.mass(stage, 1), ...
        rocket.cylinder.CoM_bor{stage, 1}, ...
        rocket.cylinder.IMoIatCoM_bor{stage, 1}] ...
        = computeCylinderProperties(rocket.cylinder, stage);
    
    % Find how far away the cylinder's center of mass is from the nosecone's
    % tip & update the distance from the nosecone's tip to the cylinder's end
    [tipToCoM, tipToOrigin] = ...
        distanceFromNoseconeTip(tipToOrigin, ...
        rocket.cylinder.CoM_bor{stage, 1}, rocket.cylinder.length(stage, 1));
    rocket.cylinder.tipToCoM_bor{stage, 1} = tipToCoM;
end

%% Total Properties
% Find the location of every component's mass center with respect to the
% nosecone's BOR coordinate system defined at the nose cone's tip with its
% x-axis pointing pointing along the nosecone's longitudinal axis towards
% its center of mass. Note that the entire dry mass rocket is aligned on
% the nosecone's centerline, so determining these positions occurs only on
% the x-axis. Also calculate each component's mass moment of inertia matrix
% with respect to the nosecone's BOR coordinate frame.
% 
% Add components to the rigid body rocket and collect the DRY:
% 1) mass
% 2) mass center
% 3) mass moment of inertia matrix,
% where the center of mass and mass moment of inertia matrix are computed
% relative to the nosecone's BOR frame.
% 
% Begin with the nosecone, which may or may not be the final staged piece.
% Regardless, count it is the final staged piece of equipment since it
% always persists (the payload is fitted into the nosecone)
rocket.dry.mass(1+rocket.stages, 1) = rocket.nosecone.mass;
rocket.dry.CoMatTip_bor{1+rocket.stages, 1} = rocket.nosecone.CoM_bor;
rocket.dry.IMoIatTip_bor{1+rocket.stages, 1} = ...
    translateIMoIatCoMToP(rocket.nosecone.mass, ...
                          rocket.nosecone.CoM_bor, ...
                          rocket.nosecone.IMoIatCoM_bor);

% Repeat the process, but now for components that aren't the nosecone
for stage = rocket.stages:-1:1
    % Get quick reference to the stage that comes after this stage (to be
    % used for obtaining properties of the next stage that have already
    % been calculated)
    nextStage = stage + 1;
    
    % Add components' masses to this stage
    newMass(1, 1) = rocket.shoulder.mass(stage, 1);
    newMass(1, 2) = rocket.cylinder.mass(stage, 1);
    % Add components' centers of mass to this stage
    newCoMs(:, 1) = rocket.shoulder.tipToCoM_bor{stage, 1};
    newCoMs(:, 2) = rocket.cylinder.tipToCoM_bor{stage, 1};
    % Add components' translated mass moments of inertia matrix to the
    % nosecone's tip
    newIMoIs(:, 1:3) = translateIMoIatCoMToP(newMass(1,1), newCoMs(:,1), rocket.shoulder.IMoIatCoM_bor{stage, 1});
	newIMoIs(:, 4:6) = translateIMoIatCoMToP(newMass(1,2), newCoMs(:,2), rocket.cylinder.IMoIatCoM_bor{stage, 1});
                          
	% Attach these components to the rocket and calculate the structural
	% properties
    [totalMass, totalCoMatTip_bor, totalIMoIatTip_bor] = ...
        attachComponents(rocket.dry.mass(nextStage, 1), ...
                         rocket.dry.CoMatTip_bor{nextStage, 1}, ...
                         rocket.dry.IMoIatTip_bor{nextStage, 1}, ...
                         newMass, newCoMs, newIMoIs);
                     
	% Add the new structural properties to the rocket's properties
    rocket.dry.mass(stage, 1) = totalMass;
    rocket.dry.CoMatTip_bor{stage, 1} = totalCoMatTip_bor;
    rocket.dry.IMoIatTip_bor{stage, 1} = totalIMoIatTip_bor;
end

% Transform the mass moment of inertia matrices from the BOR frame to the
% body-fixed principal frame. The transformation follows the form 
%                               I_t = T*I*T'
% since the mass moment of inertia matrix, I, is actually a tensor. The
% quantity I_t represents the same tensor as I but simply expressed in a
% different frame whose transformation is described by the rotation matrix
% T.
for stage = 1:rocket.stages+1
    rocket.dry.CoMatTip_body{stage, 1} = Tbody_bor*rocket.dry.CoMatTip_bor{stage, 1};
    rocket.dry.IMoIatTip_body{stage, 1} = transformTensor(Tbody_bor, rocket.dry.IMoIatTip_bor{stage, 1});
    rocket.dry.IMoIatTip_body{stage, 1} = Tbody_bor*rocket.dry.IMoIatTip_bor{stage, 1}*Tbody_bor';
end