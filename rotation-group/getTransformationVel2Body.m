function Tbod_vel = getTransformationVel2Body(angleOfAttack, angleOfVelRoll)
% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Obtain the linear transformation matrix T converting components of a
% vector expressed in the velocity frame to components of the same vector
% expressed along the body-fixed frame.
% 
%    Inputs:
% 
%     angleOfAttack - Angle of attack of the body.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    angleOfVelRoll - Angle of "velocity roll" of the body.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    Outputs:
% 
%          Tbod_vel - The transformation (rotation) matrix that converts
%                     the components of a vector expressed along the
%                     velocity axes to components of the same vector but
%                     now expressed along the body-fixed axes.
%                     Size: 3-by-3 (matrix)
%                     Units: - (unitless)
% 

% Define the rotation from the body-fixed frame to the velocity frame using
% an intermediate frame (int) realized by a rotation about the body-fixed 3
% axis through an angle `angleOfVelRoll'. Then reach the velocity frame by
% rotating about the intermediate 2-axis through an angle `angleOfAttack'.
% This process provides the rotation matrix from the body-fixed frame to
% the velocity frame, which is not so useful as the rotation matrix from
% the velocity frame to the body-fixed frame. Therefore, provide the latter
% rotation matrix instead by manually applying the transpose operation
% using negative rotation angles to rotate from the velocity frame to the
% body-fixed frame.
Tint_vel = getTransformationR2(-angleOfAttack);
Tbod_int = getTransformationR3(-angleOfVelRoll);
% Compose the rotations
Tbod_vel = Tbod_int*Tint_vel;