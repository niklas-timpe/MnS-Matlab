%% 2 a
clear all 

s = 2; %stages
delta_t = 0.1; %increment of time
tf = [0 1]; %final time
t_k = 0; %initial time

syms x t real
lambda_value = 2;
f_sym = lambda_value * x;  
f = matlabFunction(f_sym, 'Vars', {x, t});

%Butcher tableau 
A = [1/4, 1/4 - sqrt(3)/6;
     1/4 + sqrt(3)/6, 1/4];
b = [1/2, 1/2];
c = [1/2 - sqrt(3)/6, 1/2 + sqrt(3)/6];

%define K as sym and define res and jacobian as functions
syms x_k real;
syms K [length(x_k), 2] real;
[r_sym, J_sym] = symbolic_IRK_residual(f, K, x_k, t_k, delta_t, A, s);
K_func = reshape(K, [], 1);
r_func = matlabFunction(r_sym, 'Vars', {K_func, x_k});
J_func = matlabFunction(J_sym, 'Vars', {K_func, x_k});

x_k_1 = 1;
N = length(tf(1):delta_t:tf(2));
X  = zeros(N, length(x_k));
X(1,:) = x_k_1;
r = 10;
for t_k = tf(1):delta_t:tf(2)
    K_vec = repmat(x_k, numel(K), 1);
    for i = 1:5
        r = r_func(K_vec, x_k_1);
        J = J_func(K_vec,x_k_1);
        dK = -J \ r;
        K_vec = K_vec + dK;
    end
    x_k_1 = x_k_1 + delta_t * sum(b .* K_vec, 1);
    fprintf('At time t = %.2f, x = %.4f\n', t_k, x_k_1);
end



% good for now 
function [r_vec, J_sym] = symbolic_IRK_residual(f, K, x_k, t_k, delta_t, A, s)
    r_sym = sym(size(K));

    for i = 1:s
        sum_AK = sym(zeros(size(K, 1), 1));
        for j = 1:size(K, 2)
            sum_AK = sum_AK + A(i, j) * K(:, j);
        end
        
        r_sym(:, i) = f(x_k + delta_t * sum_AK, t_k) - K(:, i);
    end
    r_vec = reshape(r_sym, [], 1);
    % Reshape r_sym and K into a vector for Jacobian calculation
    J_sym = jacobian(r_vec, reshape(K, [], 1));
end