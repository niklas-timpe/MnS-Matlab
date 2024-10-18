%% 2 a 
clear all 

% Parameters
delta_t = 0.1;  % Time step
tf = 1;         % Final time

% Symbolic variable
syms x_k real;  % Define state variable symbolically
s = 2;          % Number of stages

% Define Butcher tableau
A = [1/4, 1/4 - sqrt(3)/6;
     1/4 + sqrt(3)/6, 1/4];
b = [1/2; 1/2];
c = [1/2 - sqrt(3)/6; 1/2 + sqrt(3)/6];

% Initialization
t_k = 0; 
x_k = 1;  % Initial state
Nsteps = tf / delta_t;  % Number of steps

% Loop for each time step
for k = 1:Nsteps
    % Initialize K_vec as a symbolic variable
    syms K [2, s] real;  % K should be a symbolic matrix
    K_vec = reshape(K, [], 1);  % Flatten K for use in the calculations

    iter = true;  % Initialize iter as true for the while loop
    niter = 0;    % Initialize the iteration counter

    while iter
        % Calculate the symbolic residual and Jacobian
        [r_sym, J_sym] = symbolic_IRK_residual(f, K_vec, x_k, t_k, delta_t, A, c);
        
        % Convert symbolic residual and Jacobian to Matlab functions
        r_func = matlabFunction(r_sym, 'Vars', {K_vec});
        J_func = matlabFunction(J_sym, 'Vars', {K_vec});
        
        % Evaluate the residual and Jacobian
        r = r_func(K_vec);
        J = J_func(K_vec);

        % Update K using Newton's method
        K_vec = K_vec - J \ r;  % Newton update

        % Check for convergence
        if norm(r) < 1e-5  % Adjust this tolerance if needed
            iter = false;  % Exit the loop if converged
        else
            niter = niter + 1;  % Increment the iteration counter
        end
    end

    % Reshape K_vec back to matrix form for the update step
    K_mat = reshape(K_vec, [], s);
    x_k = x_k + delta_t * (b' * K_mat);  % Compute next state using RK formula

    % Update time
    t_k = t_k + delta_t;  % Move to the next time step

    % Output the result at each step
    fprintf('At time t = %.4f, x = %.4f\n', t_k, x_k);
end

function [r_sym, J_sym] = symbolic_IRK_residual(f, K, x_k, t_k, delta_t, A, c)
    r_sym = sym(zeros(size(K)));  % Initialize residual
    sum_AK = sym(zeros(size(K, 1), 1));  % Initialize sum for each stage

    % Loop over each stage (i)
    for i = 1:size(K, 2)  % For every stage 
        sum_AK = sym(zeros(size(K, 1), 1));  % Reset sum for each stage
        for j = 1:size(K, 2)
            sum_AK = sum_AK + A(i, j) * K(:, j);  % Accumulate contributions
        end
        
        % Compute the residual for the i-th stage
        r_sym(:, i) = K(:, i) - f(x_k + delta_t * sum_AK, t_k + c(i) * delta_t);
    end
    
    r_vec = reshape(r_sym, [], 1);  % Flatten residual for Jacobian calculation
    J_sym = jacobian(r_vec, reshape(K, [], 1));  % Compute Jacobian
end 
