function Tbody_bor = getTransformationBOR2Body
% 
% Matt Werner (m.werner@vt.edu) - Jan 19, 2021
% 
% Obtain the linear transformation matrix T converting components of a
% vector expressed in the body of revolution (BOR) frame to components of
% the same vector expressed along the body-fixed frame.
% 
%    Inputs:
% 
%                   -
% 
%    Outputs:
% 
%         Tbody_bor - The transformation (rotation) matrix that converts
%                     the components of a vector expressed along the BOR
%                     axes to components of the same vector but now
%                     expressed along the body-fixed axes.
%                     Size: 3-by-3 (matrix)
%                     Units: - (unitless)
% 

% Define the transformation from the BOR frame to the body-fixed
% (principal) frame such that the body frame has its
% 3-axis (z) pointing straight out of the nosecone's tip
% 1-axis (x) lies in the cross-sectional plane of the body of revolution
% 2-axis (y) completes the right-hand rule.
% 
% This transformation comes from the orientation of the BOR frame, which is
% defined such that its
% 1-axis (x) points exactly in the opposite direction of the principal
%   frame's 3-axis (z)
% 2-axis (y) points "up" to form the usual 2D plane with x
% 3-axis (z) completes the right-hand rule.
% 
% Because the BOR is radially symmetric, the exact orientation of x and y
% in the body-fixed principal frame is rather irrelevant. Therefore, simply
% swap coordinate labels from the BOR frame to the body-fixed principal
% frame.
Tbody_bor = [0, 0, 1
             0, 1, 0
            -1, 0, 0];