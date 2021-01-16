function plotBOR(BOR, CoM)
% 
% Matt Werner (m.werner@vt.edu) - Jan 15, 2021
% 
% Plot the profile of a body of revolution (BOR). The profile consists of
% the body's cross-section containing its longitudinal axis. The result is
% the shape that would generate the 3D BOR if it were rotated 180 degrees
% about the x-axis (the x-axis is equivalent to the longitudinal axis).
% 
%    Inputs:
% 
%               BOR - Body of revolution (BOR) structure containing
%                     information about the body's sampled positions (x) at
%                     which the provided outer and inner surfaces
%                     ('outersurf' and 'innersurf', respectively) are
%                     evaluated. Also, the material density function may be
%                     provided in which case the density function is
%                     plotted beneath the BOR's profile.
% 
% ------ Overview of content contained within the BOR ------
% 
%                 x - List of sampled positions at which the outer surface,
%                     inner surface, and material density function are
%                     evaluated.
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
% 
%         outersurf - Geometric description of the body's outer surface.
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
% 
%         innersurf - Geometric description of the body's inner surface.
%                     Size: n-by-1 (vector)
%                     Units: m (meters)
% 
%           density - Optional(!) Material density function of the body
%                     along the longitudinal coordinate (x). If provided,
%                     the density function is shown as a subplot underneath
%                     the body. The density may be scalar valued, in which
%                     case it is assumed to be constant across the entire
%                     body.
%                     Size: 1-by-1 (scalar) OR n-by-1 (vector)
%                     Units: kg/m3 (kilograms per cubic meter)
% 
% ----------------------------------------------------------
% 
%               CoM - Optional(!) Position of the body's center of mass
%                     along the x-axis from the origin.
%                     Size: 3-by-1 (vector)
%                     Units: m (meters)
% 
%    Outputs:
% 
%                   -
% 

% Check to ensure that the center of mass is scalar valued
if (nargin == 2 && numel(CoM) ~= 3)
    error("Body of revolution mass center expected to be vector-valued.")
end

% Check that an appropriate number of inputs have been entered
narginchk(1, 2)

% Get the names of the fields contained within BOR
nameOfFields = string(fieldnames(BOR));

% Specify the number of subplots within the figure (1 subplot default). If
% the density is specified, then create 2 subplots within the figure.
% Regardless, the geometric outline of the body will be first.
numPlots = 1;
densityIsProvided = false;
if (any(strcmp(nameOfFields, "density")))
    numPlots = 2;
    densityIsProvided = true;
end

% Create the subplot figure
subplot(numPlots, 1, 1)

% Unpack the position and profiles of the BOR
x = BOR.x;
outersurf = BOR.outersurf;
innersurf = BOR.innersurf;

% Remove points on the inner surface along the longitudinal axis
% (especially particular for nosecones)
innersurf(innersurf < eps) = NaN;

% Plot the body
plot(x, outersurf, 'b', x, -outersurf, 'b')
hold on
plot(x, innersurf, 'r', x, -innersurf, 'r')
axis equal, grid on
% Provide a title and label the y-axis
title("Body of Revolution Cross Section", 'interpreter', 'latex')
ylabel("$y$ (m)", 'interpreter', 'latex')

% Get the x-limits
xlims = xlim;

% Check if the center of mass is provided
if (nargin == 2)
    % Extract the x- and y-components of the center of mass
    CoMx = CoM(1,1);
    CoMy = CoM(2,1);
    
    % Find the index of the closest x to the center of mass
    [~, ind] = min(abs(x - CoMx));
    
    % Determine the distance between the upper and lower inner surface at
    % the position closest to the center of mass
    verticalDistance = 2*innersurf(ind);
    
    % Convert the vertical distance to cm
    verticalDistance_cm = convUnits(verticalDistance, "m", "cm");
    
    % Mark the center of mass
    plotCenterOfMassSymbol(CoMx, CoMy, 0.068*verticalDistance_cm)
end

% Hold off on the first plot
hold off

% Check if the density is provided
if (densityIsProvided)
    % Obtain the axes for the profile plot
    profileax = gca;
    
    % Plot the density as a function of x underneath of the body
    subplot(2, 1, 2)
    
    % Unpack the density function from the BOR profile
    density = BOR.density;
    
    % Check if the density is scalar valued
    if (isscalar(density))
        % Density is assumed constant, so recast it into vector form as the
        % constant vector with size equal to that of x
        density = density*ones(size(x));
    end
    
    % Plot the material density function
    plot(x, density)
    grid on
    
    % Provide a title and label the y-axis
    title("Material Density Function", 'interpreter', 'latex')
    ylabel("$\rho$ (kg/m$^3$)", 'interpreter', 'latex')
    
    % Change the x-limits to match those of the first plot
    xlim(xlims)
    
    % Change the plotbox aspect ratio to match the first plot
    pbaspect(gca, profileax.PlotBoxAspectRatio);
end

% Plot labels
xlabel("$x$ (m)", 'interpreter', 'latex')