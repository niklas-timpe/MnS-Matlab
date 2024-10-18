%% 2 a
clear all 

s = 2; %stages
delta_t = 0.1; %increment of time
tf = 1; %final time
%x_k = 1; %initial x_k
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
K_vec = reshape(K, [], 1);
r_func = matlabFunction(r_sym, 'Vars', {K_vec, x_k});
J_func = matlabFunction(J_sym, 'Vars', {K_vec, x_k});


K_vec = zeros(numel(K), 1);

while t_k < tf
    r = rfunc(K_vec);
    J = Jfunc(K_vec);
    K_vec = K_vec - J \ r;
    K = reshape(K_vec, length(x_k), 2);
    x_k = x_k + Delta_t * x_k + Delta_t * sum(b .* K);
    t_k = t_k + delta_t;

    fprintf('At time t = %.4f, x = %.4f\n', t_k, x_k);
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

