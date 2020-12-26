function R3 = getTransformationR3(x, angleUnit)
% 
% Matt Werner (m.werner@vt.edu) - Dec 23, 2020
% 
% Calculate the elemental rotation that transfers the components of a
% vector expressed in one frame into another frame whose basis elements are
% rotated a positive (counterclockwise) angle x about the common 3-axis of
% the two frames. The result does NOT change the vector, but rather the 
% expression of its components, as in, the vector exists and applying it to
% this rotation matrix changes the frame of reference of the vector, NOT
% the vector itself.
% 
%    Inputs:
% 
%                 x - Positive angle through which the new reference frame
%                     is rotated about the two shared 3-axes of the two
%                     (new and original) frames.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians or degrees)
% 
%         angleUnit - Optional(!) Specifies whether the supplied angle x is
%                     given in units of radians or degrees.
%                     Size: 1-by-1 (string)
%                     Units: - (N/A)
% 
%                    Permissible options are:
%                    "radians" - (Default) specifies that the given angle x
%                                is in terms of radians.
% 
%                    "degrees" - (Default) specifies that the given angle x
%                                is in terms of degrees.
% 
%    Outputs:
% 
%                R3 - Rotation matrix to convert the expression of any
%                     vector from an original frame to a new frame whose
%                     basis elements are rotated a positive angle x about
%                     the common 3-axis of the two frames.
%                     Size: 3-by-3 (matrix)
%                     Units: - (N/A)
% 

% Check that the appropriate amount of inputs were supplied
narginchk(1, 2)

% Obtain the necessary sine and cosine expressions based on the supplied
% units of x
if (nargin == 1 || angleUnit == "radians")
    sinx = sin(x);
    cosx = cos(x);
elseif (angleUnit == "degrees")
    sinx = sind(x);
    cosx = cosd(x);
else
    error("Unrecognized angle unit.")
end

% Define the 3-rotation
R3 = [cosx, sinx, 0; -sinx, cosx, 0; 0, 0, 1];