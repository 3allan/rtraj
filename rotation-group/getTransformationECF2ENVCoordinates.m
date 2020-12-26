function Tenv_ecf = getTransformationECF2ENVCoordinates(longitude, geodeticLat)
% 
% Matt Werner (m.werner@vt.edu) - Dec 11, 2020
% 
% Obtain the linear transformation matrix T converting components of
% vectors expressed in ECEF coordinates to those expressed in ENV
% coordinates. Note that the ECEF and ENV frames do NOT share the same
% origin. Thus, this transformation is simply a rotation matrix between two
% frames that share the same origin, namely, the ECEF frame and another
% frame that has the same orientation as the ENV frame but NOT the same
% origin as the ENV frame.
% 
%    Inputs:
% 
%         longitude - Longitude at which the ENV frame is defined with
%                     respect to the ECF frame. This quantity is a fixed
%                     constant for a constant location on Earth's surface,
%                     such as a station or launch site.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%       geodeticLat - Geodetic latitude at which the ENV frame is defined
%                     with respect to the ECF frame. This quantity is a
%                     fixed constant for a constant location on Earth's
%                     surface, such as a station of launch site.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    Outputs:
% 
%          Tenv_ecf - The transformation (rotation) matrix that converts
%                     the components of a vector expressed along the ECF
%                     axes to components of the same vector but now
%                     expressed along the ENV axes.
%                     Size: 3-by-3 (matrix)
%                     Units: - (unitless)
% 

% Precompute trascendental quantities
sinLon = sin(longitude);
cosLon = cos(longitude);
% 
sinLat = sin(geodeticLat);
cosLat = cos(geodeticLat);

% Define the transformation matrix to rotate the coordinates from an ECF
% frame to a frame with the same orientation as the ENV frame but shares
% the same origin as the ECF frame
Tenv_ecf = [    -sinLon,         cosLon,       0;
            -sinLat*cosLon, -sinLat*sinLon, cosLat;
             cosLat*cosLon,  cosLat*sinLon, sinLat];