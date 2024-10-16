%% 2 a







%% 3 a
clc;clear;
close all
% Define the system of ODEs as an anonymous function
% x(1) = x, x(2) = y
f = @(t, x) [ 
    x(2);                    % dx/dt = y
    (1 - x(1)^2) * x(2) - x(1)  % dy/dt = (1 - x^2)y - x
];
init_conditions = [0; -1]
tf = [0 25];
options = odeset();

[t_ode45,x_ode45] = ode45(f,tf,init_conditions,options);


x_sol = x_ode45(:, 1)
y_sol = x_ode45(:, 2)

figure;
plot(t_ode45, x_sol, '-o', 'LineWidth', 1); hold on;
plot(t_ode45, y_sol, '-*', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('ODE45');
legend('x(t)', 'y(t)', 'Location', 'best');
grid on;
hold off;
saveas(gcf,'3a.png')

%% 3b

function [t, X] = rk4(f, tspan, x0, dt)
    t = tspan(1):dt:tspan(2);
    N = length(t);  % Number of time steps
    X = zeros(length(x0), N);
    X(:, 1) = x0;

    for i = 1:N-1
        tk = t(i);
        xk = X(:, i);
    
        K1 = f(tk, xk);
        K2 = f(tk, xk+dt*0.5*K1);
        K3 = f(tk, xk + 0.5*dt*K2);
        K4 = f(tk, xk+0.5*dt*K3);
    
        X(:, i+1) = xk + (dt/6) * (K1 + 2*K2 + 2*K3 + K4);
    end
end


init_conditions = [0; -1]
tf = [0 25];

[t_rk4,x_rk4] = rk4(f,tf,init_conditions,0.1);


x_sol_rk4 = x_rk4(1, :)
y_sol_rk4 = x_rk4(2, :)


delta_t = diff(t_ode45) .* 20;

figure;
plot(t_rk4, x_sol_rk4, '-*', 'Color', '#FF1D8A', 'LineWidth', 1); hold on;
plot(t_rk4, y_sol_rk4, '-', 'Color', '#47FFBE', 'LineWidth', 1);
plot(t_ode45, x_sol, '-o', 'Color', '#7E47FF', 'LineWidth', 1);
plot(t_ode45, y_sol, '-', 'Color', '#47C8FF', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('RK4 and ODE45');
plot(t_ode45(1:end-1), delta_t, 'k', 'LineWidth', 1);
legend('x_rk4(t)', 'y_rk4(t)','x_ode45(t)','y_ode45(t)', 'ode45 dt * 20', 'Location', 'best');
grid on;
hold off;
saveas(gcf,'3a3b.png')


options = odeset('AbsTol',1e-8,'RelTol',1e-8);
[t_ode45,x_ode45] = ode45(f,tf,init_conditions,options);


x_sol = x_ode45(:, 1)
y_sol = x_ode45(:, 2)

figure;
plot(t_ode45, x_sol, '.', 'Color', '#47FF83', 'LineWidth', 1); hold on;
plot(t_ode45, y_sol, '.', 'Color', '#015EFF', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('ODE45 1e-8 tolerance');
legend('x(t)', 'y(t)', 'Location', 'best');
grid on;
hold off;
saveas(gcf,'3b_ode45_tight.png')
