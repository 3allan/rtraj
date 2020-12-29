function [angleOfAttack, angleOfVelRoll] = calculateAngleOfAttack(v, V)
% 
% Matt Werner (m.werner@vt.edu) - Dec 29, 2020
% 
% Calculate the angle of attack using the velocity vector of the body as
% expressed in the body-fixed frame.
% 
%    Inputs:
% 
%                 v - Velocity of the body expressed in the body-fixed
%                     frame. This velocity is taken relative to any
%                     noninertial coordinate system attached to the Earth,
%                     that is, stationary relative to the ECF frame. Thus,
%                     the rotational velocity of Earth is NOT accounted for
%                     in this expression for the body's velocity.
%                     Size: 3-by-1 (vector)
%                     Units: m/s (meters per second)
% 
%                 V - Speed of the body expressed in the body-fixed frame.
%                     This quantity must follow the simple relation
%                                         V = |v|,
%                     which is entrusted to the user to ensure.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%    Outputs:
% 
%     angleOfAttack - Angle of attack of the body. In general, the angle of
%                     attack is nonzero and is defined as the angle made
%                     between the "forwards" direction of the body and the
%                     direction in which the body is actually travelling.
%                     It is calculated instantaneously as the angle between
%                     the vector v and the body's 3-axis, which is
%                     represented simply as [0; 0; 1] in the body-fixed
%                     frame. As such, the angle of attack is quite easy to
%                     obtain.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    angleOfVelRoll - Angle of "velocity roll" of the body, which is simply
%                     the longitude of the body-fixed vector v expressed in
%                     spherical coordinates relative to the body-fixed
%                     frame. This quantity is NOT the roll associated with
%                     the orientation of the body; it is aerodynamically
%                     unimportant, but is important with regards to
%                     obtaining a correct expression for the rotation
%                     matrices to and from the body-fixed frame and
%                     velocity frame.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 

% Obtain the angle of attack by transforming the state velocity (ENV) to
% the body-fixed frame and using the standard equation for an angle between
% two vectors

% Projection of the velocity onto the body-fixed x-y plane
vxy = [v(1:2); 0];
% Speed of the velocity projection
Vxy = norm(vxy);
% Create a perturbation to the speed if the body is precisely stationary
% relative to the ENV frame
if (V == 0), V = eps; end
if (Vxy == 0), Vxy = eps; end
% Calculate the longitudinal angle of the local velocity direction relative
% to the body-fixed frame (measured positive CCW from the body-fixed
% x-axis). Note: this quantity is NOT the kind of roll associated with the
% body's orientation relative to the chosen reference frame
angleOfVelRoll = acos(vxy(1, 1) / Vxy);
% Apply an angle correction if required
if (v(2, 1) < 0)
    angleOfVelRoll = 2*pi - angleOfVelRoll;
end
% Calculate the angle of attack
angleOfAttack = acos(v(3, 1) / V);