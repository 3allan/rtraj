function Rba = getTransformationRef2Body(q)
% 
% Matt Werner (m.werner@vt.edu)
% 
% Express the rotation matrix from a chosen reference frame to the
% body-fixed frame attached to the rigid vehicle whose orientation is
% described by the quaternion q. Applying this rotation matrix to a vector
% expressed in the reference frame produces a new expression for the SAME 
% vector, but now with its components expressed in the body-fixed frame.
% 
%    Inputs:
% 
%                 q - The supplied quaternion providing the
%                     parameterization of the body-fixed frame's
%                     orientation relative to the chosen reference frame.
%                     The quaternion follows the convention that the first
%                     three (3) elements constitute the "vector part" of
%                     the quaternion while the last (fourth/4th) element is
%                     the "scalar part" of the quaternion.
%                     Size: 4-by-1 (quaternion)
%                     Units: - (N/A)
% 
%    Outputs:
% 
%               Rba - The rotation matrix converting the components of a
%                     vector from the reference frame (a) to the body-fixed
%                     frame (b) attached to the rigid body. Applying this
%                     rotation to a vector does NOT rotate the vector
%                     relative to a single frame, but rather expresses the
%                     same vector in two different frames.
%                     Size: 3-by-3 (matrix)
%                     Units: - (N/A)
% 

% Enforce that the quaternion be unit-normalized
if (abs(norm(q) - 1) > 1e-4)
    error("Supplied quaternion is not unit-normalized.")
end

% Explicitly distribute the quaternion components for increased clarity
% --- Quaternion vector part ---
q1 = q(1);
q2 = q(2);
q3 = q(3);
% --- Quaternion scalar part ---
q4 = q(4);

% Define the rotation matrix to rotate the components of a vector from the
% reference frame to the body-fixed frame
Rba = [q1^2 - q2^2 - q3^2 + q4^2,       2*(q1*q2 + q3*q4),           2*(q1*q3 - q2*q4);
           2*(q1*q2 - q3*q4),      -q1^2 + q2^2 - q3^2 + q4^2,       2*(q2*q3 + q1*q4);
           2*(q1*q3 + q2*q4),           2*(q2*q3 - q1*q4),      -q1^2 - q2^2 + q3^2 + q4^2];