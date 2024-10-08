% Assignment 1
clear;
%% Part 1.a)
%  Coordinate System on the side view
%
%     z
%     ^
%     |
%     |
%     |______> y
%    /
%   /
%  L
% x

syms t 
syms p1(t) [3 1] 
syms phi(t) theta(t)
syms p1_dot(t) [3 1]
syms p1_dot_dot(t) [3 1]
syms phi_dot(t) 
syms theta_dot(t)
syms phi_dot_dot(t) 
syms theta_dot_dot(t)
syms q(t)
syms q_dot(t)
syms q_dot_dot(t)
syms Length
syms p2(t) [3 1]
syms g m_2 m_1
syms u [3 1]

% All q derivates with specific variables to avoid nested diffs in the
% output
q(t) = [p1(t); theta(t); phi(t)];
q_dot(t) = [p1_dot(t); theta_dot(t); phi_dot(t)];
q_dot_dot(t) = [p1_dot_dot(t); theta_dot_dot(t); phi_dot_dot(t)];

% p2 = [ a;
%        b;
%        c ]
% q = [ phi;
%       theta ]
%
% jacobian(p2,q) = [da/dphi , da/dtheta;
%                   db/dphi , db/dtheta;
%                   dc/dphi , db/dtheta]
%
% p2_dot = jacobian(p2, t) or better to read:
p2(t) = p1(t) + [-sin(phi(t))*sin(theta(t))*Length; 
           sin(phi(t))*cos(theta(t))*Length; 
           -cos(phi(t))*Length];

p2_dot(t) = jacobian(p2(t), q(t)) * q_dot(t);

% Kinteic Energy
T1 = 0.5 * m_1 * transpose(p1_dot(t)) * p1_dot(t)
T2 = 0.5 * m_2 * p2_dot(t).' * p2_dot(t);
T = T1 + T2;

% Potential Energy
V1 = [0 0 1] * m_1 * g * p1(t);
V2 = [0 0 1] * m_2 * g * p2(t);
V = V1 + V2;

% Lagrange
L = T - V;

L_grad_q_dot = gradient(T, q_dot(t));

L_grad_q_dot_diff_t = diff(L_grad_q_dot,t);
% replace diffs with variables for better readibility when printing out the
% matrices
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(phi(t)), phi_dot(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(theta(t)), theta_dot(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(phi_dot(t)), phi_dot_dot(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(theta_dot(t)), theta_dot_dot(t));

L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p11(t)), p1_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p12(t)), p1_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p13(t)), p1_dot3(t));

L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot1(t)), p1_dot_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot2(t)), p1_dot_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot3(t)), p1_dot_dot3(t));

L_grad_q = gradient(L,q(t));

EL = L_grad_q_dot_diff_t - L_grad_q;

p1_jac_q_trans = jacobian(p1(t), q(t)).';
eqn = L_grad_q_dot_diff_t - L_grad_q == p1_jac_q_trans * u;

disp("SOLUTION 1 a)")
[M,b] = equationsToMatrix(eqn,q_dot_dot(t));
disp("1 a) M Matrix = ")
disp(M)
disp("1 a) b Matrix = ")
disp(b)

% Double check the Hessian and compare it the calculated solved Weight
% Matrix
% disp("HESSIAN T = ")
% W = hessian(T, q_dot(t));
% disp(W)

%% Part 1.b)

% Probably not all clears and redefinitions are needed here but the code
% can be streamlined if needed later
clear p1(t) p2(t)
clear p1_dot(t) p2_dot(t)
clear p1_dot_dot(t) p2_dot_dot(t)
clear q(t) q_dot(t) q_dot_dot(t)
% Definition of the variables
syms z
syms p1(t) [3 1]
syms p2(t) [3 1]
syms q(t) q_dot(t) q_dot_dot(t)
syms p1_dot(t) [3 1]
syms p1_dot_dot(t) [3 1]
syms p2_dot(t) [3 1]
syms p2_dot_dot(t) [3 1]
syms C

% Generalized coordinates and derivative
q(t) = [p1; p2];
q_dot(t) = [p1_dot(t); p2_dot(t)];
q_dot_dot(t) = [p1_dot_dot(t); p2_dot_dot(t)];

% Kinetic energies
T1 = 0.5 * m_1 * transpose(p1_dot(t)) * p1_dot(t); % Helicopter
T2 = 0.5 * m_2 * transpose(p2_dot(t)) * p2_dot(t); % Hanging mass
T = T1 + T2;

% Potential energies 
V1 = m_1 * g * [0 0 1] * p1(t); % Helicopter
V2 = m_2 * g * [0 0 1] * p2(t); % Hanging mass
V = V1 + V2;

% Constraint equations 
e = p1(t) - p2(t);
C = 0.5 * (e.'* e - Length^2);

% Lagrangian (L = T - V)
L = T - V - (z*C);

% Lagrange equations
L_grad_q_dot = gradient(L, q_dot(t)); % Gradient with respect to velocities
L_grad_q = gradient(L, q(t)); % Gradient with respect to coordinates
L_grad_q_dot_diff_t = diff(L_grad_q_dot, t); % Time derivative of the gradient

% replace diffs with variables for better readibility when printing out the
% matrices
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p11(t)), p1_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p12(t)), p1_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p13(t)), p1_dot3(t));

L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot1(t)), p1_dot_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot2(t)), p1_dot_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p1_dot3(t)), p1_dot_dot3(t));

L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p21(t)), p2_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p22(t)), p2_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p23(t)), p2_dot3(t));

L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p2_dot1(t)), p2_dot_dot1(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p2_dot2(t)), p2_dot_dot2(t));
L_grad_q_dot_diff_t = subs(L_grad_q_dot_diff_t, diff(p2_dot3(t)), p2_dot_dot3(t));


p1_jac_q_trans = jacobian(p1(t), q(t)).';
eqn = L_grad_q_dot_diff_t - L_grad_q == p1_jac_q_trans * u;

disp("SOLUTION 1 b)")
[M,b] = equationsToMatrix(eqn,q_dot_dot(t));
disp("1 b) M Matrix = ")
disp(M)
disp("1 b) b Matrix = ")
disp(b)

%% Part 2.a)

% disp(eqn)
disp("SOLUTION 2 a)")
qz_vec = [q_dot_dot(t); z];

C_grad_q = gradient(C, q(t));
W = hessian(T, q_dot(t));

M = [W C_grad_q;
    C_grad_q.' 0];

disp("2 a) M Matrix = ")
disp(M)




Q = p1_jac_q_trans * u;
b1 = Q - jacobian(W*q_dot(t), q(t)) * q_dot(t) + gradient(T, q(t)) - gradient(V, q(t));
b2 = -jacobian(jacobian(C, q(t))*q_dot(t),q(t))*q_dot(t);
b = [b1;b2];

disp("2 a) b Vector = ")
disp(b)

disp("2 a) equation = ")
eqn = M * qz_vec == b;
disp(eqn)


disp("SOLUTION 2 b)")
disp("2 b) equation = ")
eqn = qz_vec == inv(M) * b;
disp(eqn)