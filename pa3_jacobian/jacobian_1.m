%% spatial jacobian
clear
clc

% const
l0 = 13;
l1 = 14.7;
l2 = 12;
l3 = 12;
l4 = 9;

eps = 1e-2;

% axes
w1 = [0, 0, 1]';
w2 = [0, -1, 0]';
w3 = [0, 0, -1]';
w4 = [0, 0, -1]';
w5 = [0, 0, -1]';
w6 = [1, 0, 0]';

% reference points
q1 = [0, -l0, l1]';
q2 = [0, -l0, l1]';
q3 = [0, -l0, l1]';
q4 = [l2, -l0, l1]';
q5 = [l2 + l3, -l0, l1]';
q6 = [l2 + l3, -l0, l1]';

% theta
theta1 = 0;
theta2 = pi/3;
theta3 = 0;
theta4 = pi/4;
theta5 = pi/3;
theta6 = pi/12;

%% exp_w
exp_w1 = compute_exp_w(w1, theta1);
exp_w2 = compute_exp_w(w2, theta2);
exp_w3 = compute_exp_w(w3, theta3);
exp_w4 = compute_exp_w(w4, theta4);
exp_w5 = compute_exp_w(w5, theta5);
exp_w6 = compute_exp_w(w6, theta6);

%% w_
w2_ = exp_w1 * w2;
w3_ = exp_w1 * exp_w2 * w3;
w4_ = exp_w1 * exp_w2 * exp_w3 * w4;
w5_ = exp_w1 * exp_w2 * exp_w3 * exp_w4 * w5;
w6_ = exp_w1 * exp_w2 * exp_w3 * exp_w4 * exp_w5 * w6;

%% q_
q4_ = q1 + exp_w1 * exp_w2 * exp_w3 * (q4 - q1);
q5_ = q4_ + exp_w1 * exp_w2 * exp_w3 * exp_w4 * (q5 - q4);
q6_ = q5_;

%% eta_
eta1 = compute_eta(w1, q1);
eta2_ = compute_eta(w2_, q2);
eta3_ = compute_eta(w3_, q3);
eta4_ = compute_eta(w4_, q4_);
eta5_ = compute_eta(w5_, q5_);
eta6_ = compute_eta(w6_, q6_);

%% Jst_s
Jst_s = [eta1, eta2_, eta3_, eta4_, eta5_, eta6_]

[U,S,V] = svd(Jst_s)

%% gst_0
gst_0_minus = [1, 0, 0, -(l2+l3+l4);
               0, 1, 0, l0;
               0, 0, 1, -l1;
               0, 0, 0, 1];

%% exp_w_minus
exp_eta6_minus = compute_exp_eta(-w6, theta6, q6);
exp_eta5_minus = compute_exp_eta(-w5, theta5, q5);
exp_eta4_minus = compute_exp_eta(-w4, theta4, q4);
exp_eta3_minus = compute_exp_eta(-w3, theta3, q3);
exp_eta2_minus = compute_exp_eta(-w2, theta2, q2);
exp_eta1_minus = compute_exp_eta(-w1, theta1, q1);

%% eta
eta1 = compute_eta(w1, q1);
eta2 = compute_eta(w2, q2);
eta3 = compute_eta(w3, q3);
eta4 = compute_eta(w4, q4);
eta5 = compute_eta(w5, q5);
eta6 = compute_eta(w6, q6);

%% eta_plus
tmp = gst_0_minus * exp_eta6_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta6_plus = Ad_g_minus * eta6;

tmp = gst_0_minus * exp_eta6_minus * exp_eta5_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta5_plus = Ad_g_minus * eta5;

tmp = gst_0_minus * exp_eta6_minus * exp_eta5_minus * exp_eta4_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta4_plus = Ad_g_minus * eta4;

tmp = gst_0_minus * exp_eta6_minus * exp_eta5_minus * exp_eta4_minus * exp_eta3_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta3_plus = Ad_g_minus * eta3;

tmp = gst_0_minus * exp_eta6_minus * exp_eta5_minus * exp_eta4_minus * exp_eta3_minus * exp_eta2_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta2_plus = Ad_g_minus * eta2;

tmp = gst_0_minus * exp_eta6_minus * exp_eta5_minus * exp_eta4_minus * exp_eta3_minus * exp_eta2_minus * exp_eta1_minus;
R = tmp(1:3, 1:3);
p = tmp(1:3, 4);
Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];
eta1_plus = Ad_g_minus * eta1;

Jst_b = [eta1_plus eta2_plus eta3_plus eta4_plus eta5_plus eta6_plus]

[S, V, D] = svd(Jst_b)