clear all;
load('input.mat')
load('output.mat')

% Split the data into estimation and validation datasets
N = length(u);
N_est = floor(N/2); % Number of samples for estimation
N_val = N - N_est;  % Number of samples for validation

% Estimation data
u_est = u(1:N_est);
y_est = y(1:N_est);

% Validation data
u_val = u(N_est+1:end);
y_val = y(N_est+1:end);


%% 3 a)

% Initialize structures to store parameters and RMSEs
theta_values = struct();
RMSE = struct();

%%% Model 12a
% y(t) + a1*y(t-1) + a2*y(t-2) = b0*u(t) + e(t)

% Make sure every run has two past values
t_est = 3:N_est;
Y_est = y_est(t_est);
Phi_est_12a = [y_est(t_est-1), y_est(t_est-2), u_est(t_est)];

% Estimate parameters using least squares
theta_12a = get_theta(Phi_est_12a, Y_est);
theta_values.model12a = struct('a1', theta_12a(1), 'a2', theta_12a(2), 'b0', theta_12a(3));

%%% Model 12b: y(t) + a1*y(t-1) + a2*y(t-2) = b0*u(t) + b1*u(t-1) + e(t)

% Make sure every run has two past values
Phi_est_12b = [y_est(t_est-1), y_est(t_est-2), u_est(t_est), u_est(t_est-1)];

% Estimate parameters using least squares
theta_12b = get_theta(Phi_est_12b, Y_est);
theta_values.model12b = struct('a1', theta_12b(1), 'a2', theta_12b(2), 'b0', theta_12b(3), 'b1', theta_12b(4));

%%% Model 12c: y(t) + a1*y(t-1) + a2*y(t-2) + a3*y(t-3) = b1*u(t-1) + e(t)

% Make sure every run has three past values
t_est_c = 4:N_est;
Y_est_c = y_est(t_est_c);
Phi_est_12c = [y_est(t_est_c-1), y_est(t_est_c-2), y_est(t_est_c-3), u_est(t_est_c-1)];

% Estimate parameters using least squares
theta_12c = get_theta(Phi_est_12c, Y_est_c);
theta_values.model12c = struct('a1', theta_12c(1), 'a2', theta_12c(2), 'a3', theta_12c(3), 'b1', theta_12c(4));

%% 3 b)

%%% Model 12a

% Prepare data for validation
t_val = 3:N_val;
Y_val = y_val(t_val);
Phi_val_12a = [y_val(t_val-1), y_val(t_val-2), u_val(t_val)];

% Predicted 1-step ahead output
Y_pred_12a = Phi_val_12a * theta_12a;

% Get RMSE
RMSE.pred_12a = sqrt(mean((Y_val - Y_pred_12a).^2));

% Simulation with validation data
y_sim_12a = zeros(N_val, 1);
% Give it initial vlaues to start of simulation
y_sim_12a(1:2) = y_val(1:2);

for t = 3:N_val
    y_sim_12a(t) = theta_values.model12a.a1 * y_sim_12a(t-1) + theta_values.model12a.a2 * y_sim_12a(t-2) + theta_values.model12a.b0 * u_val(t);
end

% Get RMSE
RMSE.sim_12a = sqrt(mean((y_val - y_sim_12a).^2));

%%% Model 12b

% Prepare data for validation
Phi_val_12b = [y_val(t_val-1), y_val(t_val-2), u_val(t_val), u_val(t_val-1)];

% Predicted 1-step ahead output
Y_pred_12b = Phi_val_12b * theta_12b;

% Get RMSE
RMSE.pred_12b = sqrt(mean((Y_val - Y_pred_12b).^2));

% Simulation with validation data
y_sim_12b = zeros(N_val, 1);
% Give it initial vlaues to start of simulation
y_sim_12b(1:2) = y_val(1:2);

for t = 3:N_val
    y_sim_12b(t) = theta_values.model12b.a1 * y_sim_12b(t-1) + theta_values.model12b.a2 * y_sim_12b(t-2) + theta_values.model12b.b0 * u_val(t) + theta_values.model12b.b1 * u_val(t-1);
end

% Get RMSE
RMSE.sim_12b = sqrt(mean((y_val - y_sim_12b).^2));

%%% Model 12c

% Prepare data for validation
t_val_c = 4:N_val;
Y_val_c = y_val(t_val_c);
Phi_val_12c = [y_val(t_val_c-1), y_val(t_val_c-2), y_val(t_val_c-3), u_val(t_val_c-1)];

% Predicted 1-step ahead output
Y_pred_12c = Phi_val_12c * theta_12c;

% Get RMSE
RMSE.pred_12c = sqrt(mean((Y_val_c - Y_pred_12c).^2));

% Simulation with validation data
y_sim_12c = zeros(N_val, 1);
% Give it initial vlaues to start of simulation
y_sim_12c(1:3) = y_val(1:3);

for t = 4:N_val
    y_sim_12c(t) = theta_values.model12c.a1 * y_sim_12c(t-1) + theta_values.model12c.a2 * y_sim_12c(t-2) + theta_values.model12c.a3 * y_sim_12c(t-3) + theta_values.model12c.b1 * u_val(t-1);
end

% Get RMSE
RMSE.sim_12c = sqrt(mean((y_val - y_sim_12c).^2));

%% FINAL Comparison

% Print all RMSE solutions
disp(RMSE)

% Get min RMSE
[~, best_sim] = min([RMSE.sim_12a, RMSE.sim_12b, RMSE.sim_12c]);
[~, best_pred] = min([RMSE.pred_12a, RMSE.pred_12b, RMSE.pred_12c]);
% Make a list of model codes for easy print out based on the above
% retrieved index
model_code = {'12a', '12b', '12c'};

disp('Best simulator: ')
disp(model_code{best_sim})

disp('Best predictor: ')
disp(model_code{best_pred})


%% FUNCTION DEFINITIONS

function theta = get_theta(phi, y_est)
    %theta = phi \ y_est;
    theta = inv(phi' * phi) * phi' * y_est
end
