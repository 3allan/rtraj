function Aprime = transformTensor(Q, A)
% 
% Matt Werner (m.werner@vt.edu) - Jan 31, 2021
% 
% Transform the representation of a (rank 2) tensor given in a particular
% coordinate frame to a new coordinate frame. The expression of rotating
% from the original coordinate system to the new coordinate system is 
% encapsulated by the information in the provided rotation matrix.
% 
%    Inputs:
% 
%                 Q - Rotation matrix to transform the provided
%                     representation of the tensor to a new coordinate
%                     system. That is, Q represents the change of basis
%                     from the original frame to the new frame.
%                     Size: 3-by-3 (matrix)
%                     Units: - (unitless)
% 
%                 A - Representation of the tensor in a particular
%                     coordinate system.
%                     Size: 3-by-3 (matrix)
%                     Units: ?
% 
%    Outputs:
% 
%            Aprime - Representation of the tensor under the basis
%                     transformation described by the rotation matrix
%                     matrix Q. That is, this matrix and the provided A
%                     matrix both represent the same tensor but simply
%                     expressed in different coordinate frames.
% 

% Transform the representation of the (rank 2) tensor to a new basis
Aprime = Q*A*Q';