% Assignment 1
clear;
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




% Not sure if we need all here but we can get rid of the unneeded once
% later
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
syms L
syms p2(t) [3 1]
syms g m_2 m_1
syms u


q = [p1; theta; phi];
q_dot = [p1_dot; theta_dot; phi_dot];
q_dot_dot = [p1_dot_dot; theta_dot_dot; phi_dot_dot];

% T = T1 + T2
% V = V1 + V2
%p1_dot = diff(p1,t);

% size 1 1
T1 = 0.5 * m_1 * transpose(p1_dot) * p1_dot


% T2 = 0.5 * m_1 * transpose(p2_dot) * p2_dot;
% size 3 1
p2 = p1 + [sin(phi)*sin(theta)*L; 
           sin(phi)*cos(theta)*L; 
           -cos(phi)*L];

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

% size 3 1
p2_dot = jacobian(p2, q) * q_dot;
disp("p2_dot = ")
disp(p2_dot)

% size 1 1
T2 = 0.5 * m_2 * p2_dot.' * p2_dot;

disp("T1 = ")
disp(T1)
disp("T1 hessian = ")
disp(hessian(T1, q_dot))
disp("T2 = ")
disp(T2)
disp("T2 hessian = ")
disp(hessian(T2, q_dot))

T = T1 + T2;
disp("HESSIAN T = ")
W = hessian(T, q_dot);
disp(W)

% size 1 1 
V1 = [0 0 1] * m_1 * g * p1;
V2 = [0 0 1] * m_2 * g * p2;
V = V1 + V2;
disp("V = ")
disp(V)

% Whole Lagrange
L = T - V;

L_grad_q_dot = gradient(T, q_dot);
disp("L gradient q_dot = ")
disp(size(L_grad_q_dot))
disp(L_grad_q_dot)
L_grad_q_dot_diff_t = diff(L_grad_q_dot,t);
disp("L gradient q_dot derivative t = ")
disp(size(L_grad_q_dot_diff_t))
disp(L_grad_q_dot_diff_t)

L_grad_q = gradient(L,q)
disp("L gradient q = ")
disp(L_grad_q)


EL = L_grad_q_dot_diff_t - L_grad_q;
disp("EL = ")
disp(EL)