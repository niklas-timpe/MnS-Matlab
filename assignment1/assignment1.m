% Assignment 1


syms q p1 theta phi
syms t L
syms x1(t) y1(t) z1(t) x2(t) y2(t) z2(t) phi(t) theta(t)
syms p1 
syms p2
syms g m_2

syms M(q) u
syms b(q, q_dot, u)

p1 = [x1; y1; z1];
p1_dot = diff(p1, t);

q = [p1; theta; phi];


% theta_dot = diff(theta, t);
% phi_dot = diff(phi, t);
% x1_dot = diff(x1, t);
% y1_dot = diff(y1,t);
% z1_dot = diff(z1, t);
% 
% x2_dot = diff(x2, t);
% y2_dot = diff(y2, t);
% z2_dot = diff(z2, t);

q_dot = diff(q, t);

p2 = p1 + [sin(phi)*cos(theta)*L; cos(phi)*L; sin(phi)*sin(theta)*L];

% q_dot = diff(q,t);
% jacob = jacobian(p2, q);
% p2_dot = jacob * q_dot

p2_dot = jacobian(p2, t)

T2 = 0.5 * m_2 * transpose(p2_dot) * p2_dot
V2 = m_2 * g * [0 0 1] * p2

v = q_dot;
v_dot = diff(v,t);

equation = M * v_dot == b;

S = solve(equation, M)