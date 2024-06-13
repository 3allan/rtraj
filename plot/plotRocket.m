function plotRocket(rocket)
% 
% Matt Werner (m.werner@vt.edu) - Feb 13, 2021
% 
% Visualize the flight vehicle (in cross-section) using its defining
% dimensions specified by the user.
% 
%            rocket - Structure containing all parameters necessary to
%                     describe a 3D flight vehicle sufficiently.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
% 

% Call a new figure
fig = figure;

% Plot the nosecone with the tip beginning at the origin and going off to
% the positive x-axis. This initial coordinate system is the nosecone BOR
% frame
x = rocket.nosecone.profile.x;
outersurf = rocket.nosecone.profile.outersurf;
innersurf = rocket.nosecone.profile.innersurf;
plot(x, outersurf, 'b', x, -outersurf, 'b')
hold on
plot(x, innersurf, 'r', x, -innersurf, 'r')
axis equal, grid on

% Assign the distance from the next BOR's origin to the nosecone
tipToOrigin = rocket.nosecone.length;

% Plot each stage
for stage = rocket.stages:-1:1
    BOR = rocket.shoulder;
    if (isstruct(BOR.profile{stage,1}))
        miniPlotBOR
    end
    
    BOR = rocket.cylinder;
    miniPlotBOR
end

% Provide a title and label the y-axis
title("Body of Revolution Cross Section", 'interpreter', 'latex')
ylabel("$y$ (m)", 'interpreter', 'latex')

function miniPlotBOR
    x = tipToOrigin + BOR.profile{stage,1}.x;
    outersurf = BOR.profile{stage,1}.outersurf;
    innersurf = BOR.profile{stage,1}.innersurf;
    plot(x, outersurf, 'b', x, -outersurf, 'b')
    plot(x, innersurf, 'r', x, -innersurf, 'r')
    tipToOrigin = tipToOrigin + BOR.length(stage,1);
end
end