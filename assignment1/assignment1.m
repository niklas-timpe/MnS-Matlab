% Assignment 1

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
syms t phi(t) theta(t)
syms p1(t) [3 1]
syms q(t)
syms L
syms p2(t) [3 1]
syms g m_2 m_1


q = [p1; theta; phi];
q_dot = diff(q, t);


% T = T1 + T2
% V = V1 + V2
p1_dot = diff(p1,t);
T1 = 0.5 * m_1 * transpose(p1_dot) * p1_dot;

% T2 = 0.5 * m_1 * transpose(p2_dot) * p2_dot;
p2 = p1 + [sin(phi)*sin(theta)*L; sin(phi)*cos(theta)*L; cos(phi)*L];

% p2_dot = jacobian(p2, t) or better to read:
p2_dot = jacobian(p2, q) * q_dot;

T2 = 0.5 * m_1 * transpose(p2_dot) * p2_dot;

T = T1 + T2;

V1 = [0 0 1] * m_1 * g * p1;
V2 = [0 0 1] * m_2 * g * p2;
V = V1 + V2;

% Whole Lagrange
L = T - V;

L_gradient_q_dot = gradient(T, q_dot) - gradient(V, q_dot);
L_gradient_q = gradient(T, q) - gradient(V, q);
L_gradient_q_dot_DOT = diff(L_gradient_q_dot, t);

simplify(L_gradient_q_dot_DOT)

