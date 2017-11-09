%% Q2
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

gst_0_minus = [1 0 0 0;
               0 1 0 -l1-l2;
               0 0 1 0;
               0 0 0 1];
           
%% exp_eta
exp_eta3_m = compute_exp_eta(-w3, theta3, q3);
exp_eta2_m = compute_exp_eta(-w2, theta2, q2);
exp_eta1_m = compute_exp_eta(-w1, theta1, q1);


%% eta3_plus
disp('eta3_plus:')

exp_eta3_minus = [ 1, 0,  0,   0;
                   0, c3, -s3, -l1*(c3 - 1);
                   0, s3, c3,  -l1*s3;
                   0, 0,  0,   1];

eta3 = [0; 0; l1; -1; 0; 0];

tmp = gst_0_minus * exp_eta3_minus;

R = tmp(1:3, 1:3);

p = tmp(1:3, 4);

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];

eta3_plus = Ad_g_minus * eta3;

eta3_plus(3) = -l2

%% eta2_plus
disp('eta2_plus:')

exp_eta2_minus = [ 1, 0,  0,   0;
                   0, c2, -s2, 0;
                   0, s2, c2,  0;
                   0, 0,  0,   1];

eta2 = [0; 0; 0; -1; 0; 0];

tmp = gst_0_minus * exp_eta3_minus * exp_eta2_minus;

R = tmp(1:3, 1:3);

p = tmp(1:3, 4);

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];

eta2_plus = Ad_g_minus * eta2;

eta2_plus(3) = -l2-l1*c3

%% eta1_plus
disp('eta1_plus:')

exp_eta1_minus = [ c1, s1, 0, 0;
                   -s1, c1,  0, 0;
                   0,  0,   1, 0;
                   0, 0,  0, 1];

eta1 = [0 0 0 0 0 1]';

tmp = gst_0_minus * exp_eta3_minus * exp_eta2_minus * exp_eta1_minus;

R = tmp(1:3, 1:3);

p = tmp(1:3, 4);

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R];

eta1_plus = Ad_g_minus * eta1;

eta1_plus(1) = -l1*c2-l2*c23;
eta1_plus(5) = -s23;
eta1_plus(6) = c23
%J = [eta1_plus, eta2_plus, eta3_plus]
