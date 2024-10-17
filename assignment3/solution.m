%% 2 a
clear all

% Define symbolic variable
syms x

% Set the value of lambda
lambda_value = -2;
% Define the dynamics as an expression with lambda substituted
f = lambda_value * x;
% Generate a MATLAB function that only depends on x
dynamics = matlabFunction(f, 'Vars', x);

tf = [0 2];
dt = 0.1;
x0 = 1;


% RK Explicit
function [ts, X] = rkbutcher(tf, dt, f, init, A, B)
    ts = tf(1):dt:tf(2);
    N = length(ts);
    s = length(B);        % Number of stages
    n = length(init);     % Dimension of the system

    X = zeros(n, N);
    X(:,1) = init;

    % Loop over simulation steps
    for step = 1:N-1
        xk = X(:, step);

        K = zeros(n, s);

        % Loop over stages
        for i = 1:s
            sum_aK = zeros(n,1);
            % Loop over stages again but only to the point that is allowed
            % for implicit euler
            for a_explicit = 1:i-1
                sum_aK = sum_aK + K(:,a_explicit) * A(i, a_explicit);
            end

            xi = xk + dt * sum_aK;
            K(:,i) = f(xi);
        end

        % Update the solution
        X(:, step+1) = xk + dt * K * B(:);
    end
end

% Define Runge-Kutta Butcher tables
Euler.A = [0];
Euler.B = [1];
Euler.C = [0];

RK2.A = [ 0  0; 
          1/2 0];
RK2.B = [0 1];
RK2.C = [0 1/2];

RK4.A = [ 0   0  0 0;
          1/2  0  0 0;
          0   1/2 0 0;
          0   0  1 0];
RK4.B = [1/6 1/3 1/3 1/6];
RK4.C = [0 1/2 1/2 0];

% Run simulations using rkbutcher
[ts_euler, X_euler] = rkbutcher(tf, dt, dynamics, x0, Euler.A, Euler.B);
[ts_rk2, X_rk2] = rkbutcher(tf, dt, dynamics, x0, RK2.A, RK2.B);
[ts_rk4, X_rk4] = rkbutcher(tf, dt, dynamics, x0, RK4.A, RK4.B);

% Plot results
figure;
plot(ts_euler, X_euler, '-o', 'DisplayName', 'Euler', 'LineWidth', 1.5);
hold on;
plot(ts_rk2, X_rk2, '-x', 'DisplayName', 'RK2', 'LineWidth', 1.5);
plot(ts_rk4, X_rk4, '-s', 'DisplayName', 'RK4', 'LineWidth', 1.5);

% Add title and labels
title('Comparison of Euler, RK2, and RK4 Methods');
xlabel('Time');
ylabel('x(t)');
legend show;
grid on;
hold off;

%% 2 b)
% Compute the true values of the system at each time step
function [X] = trueSolution(tf, dt, lambda, x0)
    ts = tf(1):dt:tf(2);
    X = zeros(1, length(ts));
    for i = 1:length(ts)
        X(i) = x0 * exp(lambda * ts(i));
    end
end

true_vals = trueSolution([0 2], 0.1, -2, 1);

% Compute and plot the error (true value - approximated value)
error_euler = true_vals - X_euler;
error_rk2 = true_vals - X_rk2;
error_rk4 = true_vals - X_rk4;

% Plot the error
figure;
plot(ts_euler, error_euler, '-o', 'DisplayName', 'Euler Error', 'LineWidth', 1.5);
hold on;
plot(ts_rk2, error_rk2, '-x', 'DisplayName', 'RK2 Error', 'LineWidth', 1.5);
plot(ts_rk4, error_rk4, '-s', 'DisplayName', 'RK4 Error', 'LineWidth', 1.5);

% Add title and labels for the error plot
title('Error Comparison of Euler, RK2, and RK4 Methods');
xlabel('Time');
ylabel('True Value - Approximated Value');
legend show;
grid on;
hold off;



%% 2 c)







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
