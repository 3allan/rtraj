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
%                     Units: ? (SI)
% 
%    Outputs:
% 
%            rocket - The newly updated structure containing all previous
%                     fields of the provided rocket plus the newly-found
%                     total mass, center of mass, and inertia matrix.
%                     Size: ? (structure)
%                     Units: ? (SI)
% 

% Required transformation to flip around axes into the body-fixed principal
% frame
Tbody_tmp = [0, 0, 1; 0, 1, 0; -1, 0, 0];

%% Nosecone
switch upper(rocket.nosecone.series)
    case "PARABOLIC"
        % Define temporary variables
        tmp_R = rocket.nosecone.OD/2;
        tmp_L = rocket.nosecone.length;
        tmp_l = rocket.nosecone.thickness/tmp_L;
        tmp_k = rocket.nosecone.k;
        
        % Get the volume
        rocket.nosecone.calculated.volume = ...
            computeParabolicConeVolume(tmp_R, tmp_L, tmp_l, tmp_k);
        % Estimate the nosecone mass
        rocket.nosecone.calculated.mass = ...
            rocket.nosecone.calculated.volume * rocket.nosecone.density;
        
        % Calculate the distance from the nosecone's tip to its center of
        % mass
        distanceTip2CoM = computeParabolicConeCenterOfMass(tmp_L, tmp_l, tmp_k);
    case "HAACK"
        
    otherwise
        error("Nosecone series not recognized.")
end

% Assign results to the official rocket structure
rocket.nosecone.calculated.mass = tmp_M;
% Transform from the temporary frame to the principal body frame
rocket.nosecone.calculated.CoM = Tbody_tmp*tmp_CoM;
rocket.nosecone.calculated.IMoI = Tbody_tmp*tmp_IMoI*Tbody_tmp';

% Clear temporary variables from the nosecone section
clear tmp_*