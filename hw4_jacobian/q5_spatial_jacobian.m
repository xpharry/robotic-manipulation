%% Q5
clear
clc

% variable
syms c1 s1 c2 s2 c3 s3 c23 s23
syms l0 l1 l2
syms theta1 theta2 theta3

w1 = [0 0 1]';
w2 = [-1 0 0]';
w3 = [-1 0 0]';

q1 = [0; 0; 0];
q2 = [0; 0; 0];
q3 = [0; l1; 0];

%% exp_w
exp_w1 = compute_exp_w(w1, theta1);
exp_w2 = compute_exp_w(w2, theta2);
exp_w3 = compute_exp_w(w3, theta3);

%% w_
w2_ = exp_w1 * w2
exp_w1 * exp_w2
w3_ = exp_w1 * exp_w2 * w3

%% eta_
eta1 = compute_eta(w1, q1)
eta2_ = compute_eta(w2_, q2)
q3_ = q1 + exp_w1 * exp_w2 * (q3-q1)
eta3_ = compute_eta(w3_, q3_)

%% Jst_s
Jst_s = [eta1, eta2_, eta3_]

J_square = Jst_s' * Jst_s

a = 1 + (l1*c2)^2 + l1^2*c1^2*s2^2 + l1^2*s1^2*s2^2

a = 1 + l1^2

det(J_square)
