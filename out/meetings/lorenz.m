% 
% Simulate the Lorenz equations
% 

close all % Closes all figure windows
clearvars % Clears workspace of all variables

%% Inputs
% Lorenz system parameters (sigma, rho, beta)
% Classical values are:
%   sigma = 10
%     rho = 28
%    beta = 8/3
% 
% See how the system changes by varying them
sigma = 10;
beta = 8/3;
rho = 80;

% ODE solver options
tspan = [0, 100];
x0 = [-8; 4; 10];
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'MaxStep', 0.05, 'Refine', 8);

%% Dynamics
% Note: The options must come before any additional parameters are passed
% in, and any parameters passed this way must be accounted for as input
% variables for the function "odefun"
sol = ode45(@odefun, tspan, x0, options, sigma, rho, beta);


%% Display
% Basic plot (with LaTeX)
figure

plot(nan, nan, 'k.', nan, nan, 'k.', nan, nan, 'k.', nan, nan, 'k.')
hold on
plot3(sol.y(1, :), sol.y(2, :), sol.y(3, :), 'Color', [0 0.4470 0.7410])
grid on
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')
title("Modified Lorenz Equations", 'Interpreter', 'latex')
view(-132, 19)
% axis tight
plot3(sol.y(1, :), sol.y(2, :), zeros(size(sol.y(3, :))), 'Color', [0.8500 0.3250 0.0980]) % xy proj
hold off
h = legend('$\sigma = 10$', '$\beta = 8/3$', '$\rho = 80$', '$\mathbf{x}(0) = (-8, 4, 10)$', 'interpreter', 'latex', 'Location', 'southoutside', 'Orientation', 'horizontal');
title(h, "Parameters", 'interpreter', 'latex')

return % End execution after first plot

% Animated plot (doesn't do any saving - just shows)
figure(2)
for k = 1:numel(sol.x)
    % Plot trajectory from IC to current step k
    plot3(sol.y(1,1:k), sol.y(2,1:k), sol.y(3,1:k))
    % Change its view
    view(40, 25)
    % Labels
    xlabel('$x$', 'Interpreter', 'latex')
    ylabel('$y$', 'Interpreter', 'latex')
    zlabel('$z$', 'Interpreter', 'latex')
    
    % Supply additional curves on the same plot
    hold on
    % Plot a red dot at the current value of the trajectory
    plot3(sol.y(1,k), sol.y(2,k), sol.y(3,k), 'r.', 'MarkerSize', 10)
    % Plot the projection onto the x-y plane
    plot3(sol.y(1,1:k), sol.y(2,1:k), zeros(1,k))
    % Plot a red dot at the current value of the trajectory in the x-y
    % plane
    plot3(sol.y(1,k), sol.y(2,k), 0, 'r.', 'MarkerSize', 10)
    % Connect the 2 red dots with a (vertical) black dashed line
    plot3(sol.y(1,k)*[1,1], sol.y(2,k)*[1,1], [0, sol.y(3,k)], 'k--')
    
    % Show grid
    grid on
    
    % Adjust axes limits
    xlim([-30, 30])
    ylim([-30, 30])
    zlim([0, 50])
    
    % Have the next iteration remove the current plot and start over
    hold off
    
    % Stop for 5 milliseconds to allow the figure to update
    pause(0.005)
end

%% Function (dxdt = f(t, x))
% State: x(t) = [x(t); y(t); z(t)]
function f = odefun(t, x, sigma, rho, beta)
    disp("t = " + t)
    f = zeros(3, 1);
    % Dynamics
    f(1) = (sigma - rho - (1 + sigma)*x(1) + x(2).*x(3))/3 + ((1 - sigma)*(x(1).^2 - x(2).^2) + 2*(sigma + rho - x(3)).*x(1).*x(2))./(3*sqrt(x(1).^2 + x(2).^2));
    f(2) = ((rho - sigma - x(3))*x(1) - (1 + sigma)*x(2))/3 + (2*(sigma - 1)*x(1).*x(2) + (sigma + rho - x(3)).*(x(1).^2 - x(2).^2))./(3*sqrt(x(1).^2 + x(2).^2));
    f(3) = (3*(x(1).^2).*x(2) - x(2).^3)/2 - beta*x(3);
end