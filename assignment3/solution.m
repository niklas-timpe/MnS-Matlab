%% Function Definitions
function [ts, X] = erkbutcher(tf, dt, f, init, A, B, C)
    ts = tf(1):dt:tf(2);
    N = length(ts);
    s = length(B);        % Number of stages
    n = length(init);     % Dimension of the system

    X = zeros(n, N);
    X(:, 1) = init;

    % Loop over simulation steps
    for step = 1:N - 1
        xk = X(:, step);
        K = zeros(n, s);

        % Loop over stages
        for i = 1:s
            sum_aK = zeros(n, 1);
            % Loop over previous stages
            for j = 1:i - 1
                sum_aK = sum_aK + A(i, j) * K(:, j);
            end
            xi = xk + dt * sum_aK;
            K(:, i) = f(ts(step) + C(i) * dt, xi);
        end

        % Update the solution
        X(:, step + 1) = xk + dt * K * B(:);
    end
end

function X = trueSolution(tf, dt, lambda, x0)
    ts = tf(1):dt:tf(2);
    X = x0 * exp(lambda * ts);
end


%% %% PART ONE
%% 2 a)

% Define symbolic variable
syms x t
lambda_value = -2;

% Define the dynamics as an expression with lambda substituted
f = lambda_value * x;
dynamics = matlabFunction(f, 'Vars', [t x]);

tf = [0 2];
dt = 0.1;
x0 = 1;

% Define Runge-Kutta Butcher tables
Euler.A = [0];
Euler.B = [1];
Euler.C = [0];

RK2.A = [0     0; 
         1/2   0];
RK2.B = [0 1];
RK2.C = [0 1/2];

RK4.A = [0     0     0   0;
         1/2   0     0   0;
         0     1/2   0   0;
         0     0     1   0];
RK4.B = [1/6 1/3 1/3 1/6];
RK4.C = [0 1/2 1/2 1];

% Run simulations using erkbutcher
[ts_euler, X_euler] = erkbutcher(tf, dt, dynamics, x0, Euler.A, Euler.B, Euler.C);
[ts_rk2, X_rk2] = erkbutcher(tf, dt, dynamics, x0, RK2.A, RK2.B, RK2.C);
[ts_rk4, X_rk4] = erkbutcher(tf, dt, dynamics, x0, RK4.A, RK4.B, RK4.C);

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

accuracy.euler = [];
accuracy.rk2 = [];
accuracy.rk4 =[];
delta_t_values = logspace(-4, 0, 1000); % This goes from 10^-3 to 10^0 (0.001 to 1)
for delta_t = delta_t_values
    % Compute the true values of the system at each time step
    true_vals = trueSolution(tf, delta_t, lambda_value, x0);
    

    [ts_euler_temp, X_euler_temp] = erkbutcher(tf, delta_t, dynamics, x0, Euler.A, Euler.B, Euler.C);
    [ts_rk2_temp, X_rk2_temp] = erkbutcher(tf, delta_t, dynamics, x0, RK2.A, RK2.B, RK2.C);
    [ts_rk4_temp, X_rk4_temp] = erkbutcher(tf, delta_t, dynamics, x0, RK4.A, RK4.B, RK4.C);
    % Compute and plot the error (true value - approximated value)
    accuracy.euler = [accuracy.euler; norm(true_vals - X_euler_temp)];
    accuracy.rk2 = [accuracy.rk2;norm(true_vals - X_rk2_temp)];
    accuracy.rk4 = [accuracy.rk4;norm(true_vals - X_rk4_temp)];
end

% Plot the error using a logarithmic scale on both axes
figure;
loglog(delta_t_values, accuracy.euler, '-', 'DisplayName', 'Euler', 'LineWidth', 1.5);
hold on;
loglog(delta_t_values, accuracy.rk2, '-', 'DisplayName', 'RK2', 'LineWidth', 1.5);
loglog(delta_t_values, accuracy.rk4, '-', 'DisplayName', 'RK4', 'LineWidth', 1.5);

% Add title and labels for the error plot
title('Error Comparison of Euler, RK2, and RK4 Methods (Log-Log Scale, Reversed X-Axis)');
xlabel('delta t (log scale, reversed)');
ylabel('Error (log scale)');
legend show;
grid on;
hold off;

true_vals = trueSolution(tf, dt, lambda_value, x0);
% Plot everything including true solution
figure;
plot(ts_euler, true_vals, 'k--', 'DisplayName', 'True Solution', 'LineWidth', 2.0);
hold on;
plot(ts_euler, X_euler, '-o', 'DisplayName', 'Euler', 'LineWidth', 1.5);
plot(ts_rk2, X_rk2, '-x', 'DisplayName', 'RK2', 'LineWidth', 1.5);
plot(ts_rk4, X_rk4, '-s', 'DisplayName', 'RK4', 'LineWidth', 1.5);

% Add title and labels
title('Comparison of Methods with True Solution');
xlabel('Time');
ylabel('x(t)');
legend show;
grid on;
hold off;

%% 3 a)
% Define the system of ODEs as an anonymous function
% x(1) = x, x(2) = y
f_sys = @(t, x) [ 
    x(2);                    % dx/dt = y
    5*(1 - x(1)^2) * x(2) - x(1)  % dy/dt = (1 - x^2)y - x
];
init_conditions = [1; 0];
tf_sys = [0 25];
options = odeset();

[t_ode45, x_ode45] = ode45(f_sys, tf_sys, init_conditions, options);

x_sol = x_ode45(:, 1);
y_sol = x_ode45(:, 2);

figure;
plot(t_ode45, x_sol, '-o', 'LineWidth', 1); hold on;
plot(t_ode45, y_sol, '-*', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('ODE45 Solutions');
legend('x(t)', 'y(t)', 'Location', 'best');
grid on;
hold off;

%% 3 b)
% Simulation parameters
dt_sys = 0.1;

% Run simulations using erkbutcher for the system
[ts_rk4_sys_m1, X_rk4_sys_m1] = erkbutcher(tf_sys, dt_sys, f_sys, init_conditions, RK4.A, RK4.B, RK4.C);

x_sol_rk4_m1 = X_rk4_sys_m1(1, :);
y_sol_rk4_m1 = X_rk4_sys_m1(2, :);

delta_t = diff(t_ode45) .* 20;

figure;
plot(ts_rk4_sys_m1, x_sol_rk4_m1, '-*', 'Color', '#FF1D8A', 'LineWidth', 1); hold on;
plot(ts_rk4_sys_m1, y_sol_rk4_m1, '-', 'Color', '#47FFBE', 'LineWidth', 1);
plot(t_ode45, x_sol, '-o', 'Color', '#7E47FF', 'LineWidth', 1);
plot(t_ode45, y_sol, '-', 'Color', '#47C8FF', 'LineWidth', 1);
plot(t_ode45(1:end-1), delta_t, 'k', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('RK4 and ODE45 Solutions');
legend('x_{RK4}(t)', 'y_{RK4}(t)', 'x_{ODE45}(t)', 'y_{ODE45}(t)', 'ODE45 dt * 20', 'Location', 'best');
grid on;
hold off;



dt_sys = 0.01;

[ts_rk4_sys, X_rk4_sys] = erkbutcher(tf_sys, dt_sys, f_sys, init_conditions, RK4.A, RK4.B, RK4.C);

x_sol_rk4 = X_rk4_sys(1, :);
y_sol_rk4 = X_rk4_sys(2, :);

% ODE45 with tighter tolerances
options_tight = odeset('AbsTol',1e-8,'RelTol',1e-8);
[t_ode45_tight, x_ode45_tight] = ode45(f_sys, tf_sys, init_conditions, options_tight);

x_sol_tight = x_ode45_tight(:, 1);
y_sol_tight = x_ode45_tight(:, 2);

figure;
plot(ts_rk4_sys, x_sol_rk4, '-o', 'Color', '#FF1D8A', 'LineWidth', 1); hold on;
plot(ts_rk4_sys, y_sol_rk4, '-o', 'Color', '#47FFBE', 'LineWidth', 1);
plot(t_ode45_tight, x_sol_tight, '-*', 'Color', '#47FF83', 'LineWidth', 1);
plot(t_ode45_tight, y_sol_tight, '-*', 'Color', '#015EFF', 'LineWidth', 1);
xlabel('Time t');
ylabel('Solutions');
title('ODE45 with Tight Tolerances vs RK4 with 0.01 dt');
legend('x_{RK4}(t)', 'y_{RK4}(t)', 'x_{ODE45}(t)', 'y_{ODE45}(t)', 'Location', 'best');
grid on;
hold off;

%% %% PART TWO
%% 4 a)
% Define symbolic variables and functions for IRK4
n = 1;
s = 2;

% Create Input Function
x_sym = sym('x', [n, 1]);
syms t_sym

% 4 b)
f_sym = -2 * x_sym;

% Generate InputFunction.m file
matlabFunction(f_sym, 'file', 'InputFunction', 'vars', {t_sym, x_sym});

% Define the Butcher Tableau for IRK4
A_IRK4 = [1/4, 1/4 - sqrt(3)/6; 
          1/4 + sqrt(3)/6, 1/4];
c_IRK4 = [1/2 - sqrt(3)/6; 
          1/2 + sqrt(3)/6];
b_IRK4 = [1/2, 1/2];

% Defining Residual Function
K_mat_sym = sym('K', [n, s]);
r_mat_sym = sym('r', [n, s]);
syms t dt real

for i = 1:s
    r_mat_sym(:, i) = InputFunction(t + c_IRK4(i) * dt, x_sym + dt * K_mat_sym * A_IRK4(i, :)') - K_mat_sym(:, i);
end

% Reshape K and r to vectors to enable computation of Jacobian
r_vec = reshape(r_mat_sym, [n * s, 1]);
K_vec = reshape(K_mat_sym, [n * s, 1]);

dr = jacobian(r_vec, K_vec);

% Generate IRK4.m file
matlabFunction(r_vec, dr, 'file', 'IRK4', 'vars', {t, x_sym, K_vec, dt});

% Simulation parameters:
tf_IRK4 = [0 1];
dt_IRK4 = 0.1;
x0_IRK4 = [1];

% Simulate using IRK4
Nsteps_IRK4 = (tf_IRK4(2) - tf_IRK4(1)) / dt_IRK4;
t_IRK4 = tf_IRK4(1):dt_IRK4:tf_IRK4(2);
x_IRK4 = zeros(n, length(t_IRK4));
x_IRK4(:, 1) = x0_IRK4;

% Loop for the IRK4
for k = 1:Nsteps_IRK4
    % Newton iteration
    iter = true;
    niter = 0;
    % Initialize K_i = x_k
    K_vec = repmat(x_IRK4(:, k), s, 1);

    while iter
        [r_val, dr_val] = IRK4(t_IRK4(k), x_IRK4(:, k), K_vec, dt_IRK4);
        delta_K = dr_val \ r_val;
        K_vec = K_vec - delta_K;

        if norm(r_val) < 1e-5 || niter > 100
            iter = false;
        else
            niter = niter + 1;
        end
    end
    % Reshape K to matrix for update step
    K_mat = reshape(K_vec, [n, s]);
    x_IRK4(:, k + 1) = x_IRK4(:, k) + dt_IRK4 * K_mat * b_IRK4';
end



figure;
plot(t_IRK4, x_IRK4, '-o', 'DisplayName', 'IRK4', 'LineWidth', 1.5);hold on;
plot(ts_rk4, X_rk4, '*','DisplayName', 'RK4', 'LineWidth', 1.5)
title('IRK4 Method vs RK4 from 2 Solution');
xlabel('Time');
ylabel('x(t)');
legend show;
grid on;
hold off;




%% 4 c)
% Create Input Function
% Define symbolic variables and functions for IRK4
n = 2;
s = 2;

x_sym = sym('x', [n, 1]);
syms t_sym
vdp = [ 
    x_sym(2);                                 % dx/dt = y
    5*(1 - x_sym(1)^2) * x_sym(2) - x_sym(1)  % dy/dt = (1 - x^2)y - x
];


% Generate InputFunction.m file
matlabFunction(vdp, 'file', 'VanDerPol', 'vars', {t_sym, x_sym});

% Define the Butcher Tableau for IRK4
A_IRK4 = [1/4, 1/4 - sqrt(3)/6; 
          1/4 + sqrt(3)/6, 1/4];
c_IRK4 = [1/2 - sqrt(3)/6; 
          1/2 + sqrt(3)/6];
b_IRK4 = [1/2, 1/2];

% Defining Residual Function
K_mat_sym = sym('K', [n, s]);
r_mat_sym = sym('r', [n, s]);
syms t dt real

for i = 1:s
    r_mat_sym(:, i) = VanDerPol(t + c_IRK4(i) * dt, x_sym + dt * K_mat_sym * A_IRK4(i, :)') - K_mat_sym(:, i);
end

% Reshape K and r to vectors to enable computation of Jacobian
r_vec = reshape(r_mat_sym, [n * s, 1]);
K_vec = reshape(K_mat_sym, [n * s, 1]);

dr = jacobian(r_vec, K_vec);

% Generate IRK4.m file
matlabFunction(r_vec, dr, 'file', 'IRK4', 'vars', {t, x_sym, K_vec, dt});

% Simulation parameters:
tf_IRK4 = [0 25];
dt_IRK4 = 0.01;
x0_IRK4 = [1;0];

% Simulate using IRK4
Nsteps_IRK4 = (tf_IRK4(2) - tf_IRK4(1)) / dt_IRK4;
t_IRK4 = tf_IRK4(1):dt_IRK4:tf_IRK4(2);
x_IRK4 = zeros(n, length(t_IRK4));
x_IRK4(:, 1) = x0_IRK4;

% Loop for the IRK4
for k = 1:Nsteps_IRK4
    % Newton iteration
    iter = true;
    niter = 0;
    % Initialize K_i = x_k
    K_vec = repmat(x_IRK4(:, k), s, 1);

    while iter
        [r_val, dr_val] = IRK4(t_IRK4(k), x_IRK4(:, k), K_vec, dt_IRK4);
        delta_K = dr_val \ r_val;
        K_vec = K_vec - delta_K;

        if norm(r_val) < 1e-5 || niter > 100
            iter = false;
        else
            niter = niter + 1;
        end
    end
    % Reshape K to matrix for update step
    K_mat = reshape(K_vec, [n, s]);
    x_IRK4(:, k + 1) = x_IRK4(:, k) + dt_IRK4 * K_mat * b_IRK4';
end



figure;
plot(ts_rk4_sys, x_sol_rk4, '-o', 'Color', '#FF1D8A', 'LineWidth', 1); hold on;
plot(ts_rk4_sys, y_sol_rk4, '-o', 'Color', '#47FFBE', 'LineWidth', 1);
plot(t_IRK4, x_IRK4(1,:), '-*', 'Color', '#7E47FF', 'LineWidth', 1);
plot(t_IRK4, x_IRK4(2,:), '-*', 'Color', '#47C8FF', 'LineWidth', 1);
legend('x_{RK4}(t)', 'y_{RK}(t)', 'x_{IK4}(t)', 'y_{IRK4}(t)', 'Location', 'best');
title('IRK4 VDP vs Butcher');
xlabel('Time');
ylabel('x(t)');
legend show;
grid on;
hold off;
