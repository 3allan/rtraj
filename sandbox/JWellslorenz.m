%
% Simulate the Lorenz equations
%
clear

sigma = 10;
rho = 80;
beta = 8/3;

tspan = [0, 100];
x0 = [-8; 4; 10];

options = odeset('RelTol', 1e-10, 'Refine', 6);
%tried lower tolerances, no difference in time and looked worse

sol = ode45(@odefun, tspan, x0, options, sigma, rho, beta);

figure(1)
plot3(sol.y(1, :), sol.y(2, :), sol.y(3, :))
grid on
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
zlabel('$z$', 'Interpreter', 'latex')
view(-132, 19)


figure(2)

for k = 1:numel(sol.x)
    %tried to make plot3 only plot the new point, but that went very wrong
    plot3(sol.y(1, 1:k), sol.y(2, 1:k), sol.y(3, 1:k))
    %x,y,z labels by far take the most amount of time, made it so they only
    %show up at the end
    if k == numel(sol.x)
        xlabel('$x$', 'Interpreter', 'latex')
        ylabel('$y$', 'Interpreter', 'latex')
        zlabel('$z$', 'Interpreter', 'latex')
    end
    view(-132, 19)
    
    hold on
    plot3(sol.y(1, k), sol.y(2, k), sol.y(3, k), 'r.', 'MarkerSize', 10)
    
    plot3(sol.y(1, 1:k), sol.y(2, 1:k), zeros(1,k))
    plot3(sol.y(1, k), sol.y(2, k), 0, 'r.', 'MarkerSize', 10)
    plot3(sol.y(1, k)*[1,1], sol.y(2, k)*[1,1], [0, sol.y(3, k)], 'k--')
    
    grid on
    
    xlim([-20 20])
    ylim([-10 10])
    zlim([0 150])
    
    pause(0.001)
    hold off
end

%% Function (dxdt = f(t,x))
% State x(t) = [x(t); y(t); z(t)]
function f = odefun(~,x, sigma, rho, beta)
    f = zeros(3, 1);
    %Dynamics
    %tried saving time by assigning variable to limit number of
    %calculations, but that took more time
    f(1) = (sigma-rho-(1+sigma)*x(1)+x(2)*x(3))/3+((1-sigma)*(x(1)^2-x(2)^2)+2*(sigma+rho-x(3))*x(1)*x(2))/(3*sqrt(x(1)^2+x(2)^2));
    f(2) = ((rho-sigma-x(3))*x(1)-(sigma+1)*x(2))/3+(2*(sigma-1)*x(1)*x(2)+(sigma+rho-x(3))*(x(1)^2-x(2)^2))/(3*sqrt(x(1)^2+x(2)^2));
    f(3) = (3*x(1)^2*x(2)-x(2)^3)/2 - beta*x(3);
end
