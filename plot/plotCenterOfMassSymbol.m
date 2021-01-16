function plotCenterOfMassSymbol(X, Y, radius, ax)
% 
% Matt Werner (m.werner@vt.edu) - Jan 16, 2021
% 
% Plot the alternating black and white quarter sectors of a circle commonly
% used to mark the center of mass of an object. This symbol is 2D.
% 
%    Inputs:
% 
%                 X - Position along the x-axis that the symbol's center is
%                     located.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 Y - Position along the y-axis that the symbol's center is
%                     located.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%            radius - Radius of the center of mass symbol.
%                     Size: 1-by-1 (scalar)
%                     Units: cm (centimeters)
% 
%                ax - Optional(!) Axes of the figure to which the symbol
%                     should be added. If left unspecified, then the
%                     default behavior is to add the symbol to gca.
%                     Size: 1-by-1 (axes)
%                     Units: N/A
% 

% Check inputs
narginchk(3, 4)

% Add the symbol to gca if a set of axes aren't specified
if (nargin == 3)
    ax = gca;
end

% Obtain the radius in meters
convertedRadius = convUnits(radius, "cm", "m");

% Parameterize the main black circle that exists in the background behind
% the two white sectors that will be filled overtop of it
t = linspace(0, 1, 500);
[x, y] = arcCoords;

% Plot the black background circle
fill(ax, X+x, Y+y, 'k')
hold on

% Update the parameterization for the upper-right quarter sector of the
% circle while leaving some room for a black border around the white
% sectors
t = t/4;
[x, y] = arcCoords;
x = 0.98*[0, x, 0];
y = 0.98*[0, y, 0];

% Plot the bottom-right white sector
fill(ax, X+x, Y-y, 'w')
% Plot the top-left white sector
fill(ax, X-x, Y+y, 'w')
hold off

function [x, y] = arcCoords
% Obtain the x and y coordinates of an arc of the circle parameterized by
% the argument t and radius 'convertedRadius', whose units are meters.
twopit = 2*pi*t;
x = convertedRadius*cos(twopit);
y = convertedRadius*sin(twopit);
end
end