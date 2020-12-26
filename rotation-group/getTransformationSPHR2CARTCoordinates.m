function Tcart_sphr = getTransformationSPHR2CARTCoordinates(longitude, colatitude)
% 
% Matt Werner (m.werner@vt.edu) - Dec 11, 2020
% 
% Obtain the linear transformation matrix T converting components of a
% vector expressed along a spherical frame (defined at a longitude and
% geocentric latitude with respect to a Cartesian frame) to components of
% the same vector expressed along the Cartesian frame. Note that,
% physically, these two frames do NOT share the same origin in general.
% Thus, this transformation is simply a rotation matrix.
% 
%    Inputs:
% 
%         longitude - Longitude at which the spherical frame is defined
%                     with respect to the Cartesian frame.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%        colatitude - Spherical colatitude at which the spherical frame is
%                     defined with respect to the Cartesian frame. Note
%                     that this quantity is NOT the latitude. The two are
%                     related by
%                                latitude + colatitude = pi/2
%                     such that
%                                sin(colatitude) = cos(latitude)
%                                cos(colatitude) = sin(latitude).
%                     Spherical coordinates are very often reported with
%                     respect to the colatitude as opposed to latitude.
%                     Size: 1-by-1 (scalar)
%                     Units: - (radians)
% 
%    Outputs:
% 
%        Tcart_sphr - The transformation (rotation) matrix that converts
%                     the components of a vector expressed along the
%                     spherical frame defined some distance away from the
%                     Cartesian frame's origin at a longitude and
%                     colatitude to components of the same vector, but
%                     expressed along the Cartesian frame's axes.
%                     Size: 3-by-3 (matrix)
%                     Units: - (unitless)

% Precompute trascendental quantities
sinLon = sin(longitude);
cosLon = cos(longitude);
% 
sinColat = sin(colatitude);
cosColat = cos(colatitude);

% Define the transformation matrix to rotate the coordinates from an ECF
% frame to a frame with the same orientation as the spherical frame formed
% at a geocentric latitude and longitude (arbitrary distance away from the
% ECF origin)
Tcart_sphr = [sinColat*cosLon, cosColat*cosLon, -sinLon;
              sinColat*sinLon, cosColat*sinLon, +cosLon;
                  cosColat,        -sinColat,       0  ];