% 
% Simulate the Lorenz equations
% 

close all % Closes all figure windows
clearvars % Clears workspace of all variables
clc % Clears command window

%% Inputs
% Lorenz system parameters (sigma, rho, beta)
% Classical values are:
%   sigma = 10
%     rho = 28
%    beta = 8/3
% 
%Toggle use of animation graph (takes a really long time)
animating = true;

% See how the system changes by varying them
sigma = 10;
rho = 80;
beta = 8/3;

% ODE solver options
tspan = [0, 100];
x0 = [-8; 4; 10];
options = odeset('RelTol', 1e-13,'MaxStep', 0.05, 'Refine', 4);

%% Dynamics
% Note: The options must come before any additional parameters are passed
% in, and any parameters passed this way must be accounted for as input
% variables for the function "odefun"
sol = ode45(@odefun, tspan, x0, options, sigma, rho, beta);

%% Basic Display
% Basic plot (with LaTeX)
figure(1)
plot3(sol.y(1, :), sol.y(2, :), sol.y(3, :))
grid on
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')
view(-132, 19)
hold off

% return % End execution after first plot

%% Animated plot (doesn't do any saving - just shows)
if animating
    figure(2)
    % Animated line for trajectory in 3d (blue)
    traj3d = animatedline('Color',[0 0.4470 0.7410]);
    % Animated line for projection of trajectory into x-y plane (orange)
    traj2d = animatedline('Color',[0.9290 0.6940 0.1250]);
    % Animated line for red dot on current point in 3d trajectory
    dot3d = animatedline('Color','r','Marker','.','MarkerSize',10,'MaximumNumPoints',1);
    % Animated line for red dot on current point in 2d projection
    dot2d = animatedline('Color','r','Marker','.','MarkerSize',10,'MaximumNumPoints',1);
    % Animated line between two red dots (black dashed line)
    vertLine = animatedline('Color','k','LineStyle','--','MaximumNumPoints',2);

    % PLOT FORMATTING
    % Change its view
    view(-132, 19)
    % Show grid
    grid on
    % Labels
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    zlabel('$z$', 'Interpreter', 'latex')
    % Adjust axes limits
    xlim([-20, 20])
    ylim([-10, 10])
    zlim([0, 150])

    %PLOT ANIMATION
    for k = 1:numel(sol.x)
        %Add current point to 3d trajectory
        addpoints(traj3d,sol.y(1,k), sol.y(2,k), sol.y(3,k))
        %Add current x,y coordinates to 2d projection of trajectory
        addpoints(traj2d,sol.y(1,k), sol.y(2,k), 0)
        %Move red dot to current point in 3d
        addpoints(dot3d,sol.y(1,k), sol.y(2,k), sol.y(3,k))
        %Move red dot to current x,y coordinates in 2d projection
        addpoints(dot2d,sol.y(1,k), sol.y(2,k),0)
        %Move vertical line to connect current red dots
        addpoints(vertLine,sol.y(1,k)*[1 1], sol.y(2,k)*[1 1],sol.y(3,k)*[0 1]);

        %Brief pause for animation to be updated
        pause(0.00001)

    end
end
%% Dynamic Reconstruction
% Pass solution into function f to obtain time derivatives of x, y, and z
xdot = odefun(tspan,sol.y,sigma,rho,beta);
%table(sol.x',xdot(1,:)',xdot(2,:)',xdot(3,:)','VariableNames',["t","xdot","ydot","zdot"])

%% Function (dxdt = f(t, x))
% State: x(t) = [x(t); y(t); z(t)]
function f = odefun(~, x, sigma, rho, beta)
    f = zeros(3, length(x(1,:)));
    % Dynamics
    %f(1) = sigma*(x(2) - x(1));
    %f(2) = x(1)*(rho - x(3)) - x(2);
    %f(3) = x(1)*x(2) - beta*x(3);
    f(1,:) = 1/3 .* (sigma - rho - (1 + sigma).*x(1,:) + x(2,:).*x(3,:)) + ...
        ((1-sigma).*(x(1,:).^2 - x(2,:).^2) + 2.*(sigma + rho - x(3,:)).*x(1,:).*x(2,:))./(3.*sqrt(x(1,:).^2 + x(2,:).^2));
    f(2,:) = 1/3 .* ((rho - sigma - x(3,:)).*x(1,:) - (sigma + 1).*x(2,:)) + ...
        (2.*(sigma - 1).*x(1,:).*x(2,:) + (sigma + rho - x(3,:)) .* (x(1,:).^2 - x(2,:).^2))./(3.*sqrt(x(1,:).^2 + x(2,:).^2));
    f(3,:) = 1/2 .* (3.*x(1,:).^2.*x(2,:) - x(2,:).^3) - beta.*x(3,:);
end