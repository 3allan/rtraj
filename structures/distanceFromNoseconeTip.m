function [tipToCoM, tipToEnd] = ...
    distanceFromNoseconeTip(tipToOrigin, originToCoM, componentLength)
% 
% Matt Werner (m.werner@vt.edu) - Jan 31, 2021
% 
% Calculate how far away a component's center of mass as projected onto the
% rocket's centerline is from the nosecone's tip. Further, calculate the
% complete distance from the nosecone's tip to the full extent of the
% component.
% 
%    Inputs:
% 
%       TipToOrigin - Distance from the nosecone's tip to the origin of
%                     this component's BOR frame along the vehicle's
%                     centerline.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%        Origin2CoM - Distance from the origin of this component's BOR
%                     frame to this component's center of mass as along the
%                     vehicle's centerline.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%   componentLength - Distance that this component spans along the rocket's
%                     centerline.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%    Outputs:
% 
%          TipToCoM - Distance from the nosecone's tip to the center of
%                     mass of this component.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%          TipToEnd - Distance from the nsoecone's tip to the 'end' of this
%                     component.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 

% Calculate the distance from the nosecone's tip to the component's center
% of mass
tipToCoM = tipToOrigin + originToCoM;

% Calculate the distance from the nosecone's tip to the extent of the
% component
tipToEnd = [tipToOrigin(1,1) + componentLength; 0 ;0];