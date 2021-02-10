function rocket = computeWetProperties(rocket)
% 
% Matt Werner (m.werner@vt.edu) - Feb 8, 2021
% 
% Determine the wet mass, wet center of mass, and wet mass moment of
% inertia matrix for the flight vehicle. Here, 'wet' means fully loaded
% such that, at each stage, the properties are calculated before any
% fuel/propellant has been utilized.
% 
% These values are calculated first in the nosecone tip's BOR frame and
% then transformed to be expressed in the nosecone tip's body (principal)
% frame.
% 
% The motors are modelled as BATES grain cylinders for the purpose of
% determining their center of mass and mass moment of inertia. The center
% of mass remains fixed in the exact middle of the motor until it vanishes
% when the fuel/motor is fully gone. If a burn depth profile is not given,
% then the motor is assumed to be a solid (ie not hollow) cylinder varying
% in web thickness from 0 to ID/2, where ID is the casing's inner diameter.
% 
%    Inputs:
% 
%            rocket - The input details of the rocket in standard SI units
%                     concerning the fuel tanks/motors. These parameters
%                     are used in determining the rocket's wet structural 
%                     properties (mass, center of mass, and mass moment of
%                     inertia matrix).
%                     Size: ? (structure)
%                     Units: SI
% 
%    Outputs:
% 
%            rocket - The newly updated rocket containing information
%                     regarding mass, center of mass, and mass moment of
%                     inertia matrix of the wet vehicle. These quantities
%                     are taken with respect to the nosecone tip's BOR and
%                     body (principal) frames. This vehicle is fully
%                     "assembled" such that its structure is rigidly
%                     attached together including the fuel tanks/motors.
%                     This vehicle can be plotted to show its geometry.
%                     Size: ? (structure)
%                     Units: SI
% 

% Required transformation to flip the BOR coordinate frame into the
% body-fixed principal frame
Tbody_bor = getTransformationBOR2Body;

%% Motor Properties & All Wet Masses
% "Calculate" the masses first to ensure that all motors are ready for
% calculating the center of mass and mass moments of inertia
for stage = 1:rocket.stages
    % ===================== Ensure Proper Sizing ==========================
    % Check that the motor fits into the casing lengthwise
    if (rocket.motor.length(stage,1) > rocket.cylinder.length(stage,1))
        error("Motor does not fit into casing (stage %1.0f)", stage)
    end
    
    % ================ Obtain Properties in its BOR Frame =================
    % Evaluate the structural properties of the motor
    [rocket.motor.profile{stage,1}, ...
        rocket.motor.mass(stage,1), ...
        rocket.motor.CoM_bor{stage,1}, ...
        rocket.motor.IMoIatCoM_bor{stage,1}, ...
        rocket.motor.dIMoIatCoMdt_bor{stage,1}] = ...
            computeMotorProperties(rocket.motor, stage);
    
    % ============================= Mass ==================================
    % Calculate the rocket's wet mass (all remaining dry mass + motors)
    rocket.wet.mass(stage, 1) = rocket.dry.mass(stage, 1) + ...
                      sum(rocket.motor.mass(stage:end, 1));
end

%% Wet Center of Mass
% Use the dry properties to calculate the wet properties
for stage = 1:rocket.stages
    % ========================= Center of Mass (BOR) ============================
    % Collect all of the motors' centers of mass into a matrix prior to the
    % calculation and also collect the motors' masses.
    for motor = rocket.stages:-1:stage
        % Intermediate variables for the calculation
        tmp_cylinder_tipToCoMs(:, motor) = rocket.cylinder.tipToCoM_bor{motor, 1};
        tmp_cylinder_CoMs(:, motor) = rocket.cylinder.CoM_bor{motor, 1};
        tmp_cylinder_lengths(:, motor) = [rocket.cylinder.length(motor, 1); 0; 0];
        tmp_motor_CoMs(:, motor) = rocket.motor.CoM_bor{motor, 1};
        % Prepare the masses as well
        tmp_motor_masses(1, motor) = rocket.motor.mass(motor, 1);
    end
    % Note: Leave room here to account for the nozzles ~~~~~~~~~~~~~~~~~~~~~~~~~~
    tmp_motor_tipToCoMs = tmp_cylinder_tipToCoMs - tmp_cylinder_CoMs + tmp_cylinder_lengths - ...
        tmp_motor_CoMs;
    % Assign the distance from nosetip to center of mass for this motor
    rocket.motor.tipToCoM_bor{stage, 1} = tmp_motor_tipToCoMs(:, stage);
    % Calculate the location of the motor's center of mass with respect to
	% the nosecone tip's BOR frame
    rocket.wet.CoMatTip_bor{stage,1} = ...
        (rocket.dry.CoMatTip_bor{stage,1} * rocket.dry.mass(stage,1) + ...
         sum(tmp_motor_tipToCoMs.*tmp_motor_masses, 2))/rocket.wet.mass(stage,1);
     clear tmp_*
end

%% Translate All Moments of Inertia to Nosecone Tip
% Translate all of the motors' mass moments of inertia over to the nosecone
% tip BOR frame
for stage = 1:rocket.stages
    % ================ Mass Moment of Inertia Matrix (BOR) ================
    % Translate the motor's mass moment of inertia to the nosecone tip in
    % the BOR frame
    rocket.motor.IMoIatTip_bor{stage, 1} = ...
        translateIMoIatCoMToP(rocket.motor.mass(stage, 1), ...
        rocket.motor.tipToCoM_bor{stage,1}, ...
        rocket.motor.IMoIatCoM_bor{stage,1});
end

%% Place Motors Inside of the Casings
% Calculate the total wet moment of inertia matrix (BOR frame) for each
% fully loaded stage in the rocket
for stage = 1:rocket.stages
    % Initial wet mass moment of inertia
    for motor = stage:rocket.stages
         rocket.wet.IMoIatTip_bor{stage, 1} = ...
             rocket.dry.IMoIatTip_bor{stage, 1} + ...
             rocket.motor.IMoIatTip_bor{stage, 1};
    end
end

%% Allocate the Nosecone as the Final Possible Stage
% Include the nosecone as the final stage in case the mission requires that
% it be separated from the final sustainer
rocket.wet.mass(1+rocket.stages,1) = rocket.dry.mass(1+rocket.stages,1);
rocket.wet.CoMatTip_bor{1+rocket.stages, 1} = rocket.dry.CoMatTip_bor{1+rocket.stages, 1};
rocket.wet.IMoIatTip_bor{1+rocket.stages, 1} = rocket.dry.IMoIatTip_bor{1+rocket.stages, 1};

%% Basis Transformation of the Centers of Mass and Inertia Matrices
% Transform the mass moment of inertia matrices from the BOR frame to the
% body-fixed principal frame. The transformation follows the form 
%                           I_t = T*I*T'
% since the mass moment of inertia matrix, I, is actually a tensor. The
% quantity I_t represents the same tensor as I but simply expressed in a
% different frame whose transformation is described by the rotation matrix
% T.
for stage = 1:rocket.stages
    % Transform motor BOR properties into the body (principal) frame
    % Note: The motor is the only rocket component that gets its center of
    % mass and inertia matrix transformed into the body-fixed principal
    % frame since it is the only component that continuously changes in
    % time (thus affecting the actual center of mass and inertia matrix of
    % the rocket throughout the entire burn time)
    rocket.motor.CoM_body{stage,1} = Tbody_bor*rocket.motor.CoM_bor{stage,1};
    rocket.motor.IMoIatCoM_body{stage,1} = transformTensor(Tbody_bor, rocket.motor.IMoIatCoM_bor{stage,1});
    rocket.motor.dIMoIatCoMdt_body{stage,1} = transformTensor(Tbody_bor, rocket.motor.dIMoIatCoMdt_bor{stage,1});
    rocket.motor.tipToCoM_body{stage,1} = Tbody_bor*rocket.motor.tipToCoM_bor{stage,1};
    % Transform wet BOR properties into the body (principal) frame
    rocket.wet.CoMatTip_body{stage,1} = Tbody_bor*rocket.wet.CoMatTip_bor{stage,1};
    rocket.wet.IMoIatTip_body{stage, 1} = transformTensor(Tbody_bor, rocket.wet.IMoIatTip_bor{stage,1});
end
