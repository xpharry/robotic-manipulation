%% Q2
clear
clc

% variable
syms c1 s1 c2 s2 c3 s3
syms l1

q1 = [0, 0, 0];
q2 = [0, 0, 0];
q3 = [0, 0, 0];

gst_0_minus = [1 0 0 0;
               0 1 0 -l1;
               0 0 1 0;
               0 0 0 1];
           
%% eta3_plus

exp_eta3_minus = [c3 0 -s3 0;
                  0 1 0 0;
                  s3 0 c3 0;
                  0 0 0 1];

eta3 = [0 0 0 0 1 0]';

tmp = gst_0_minus * exp_eta3_minus;

R = tmp(1:3, 1:3);

p = tmp(1:3, 4);

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R]

eta3_plus = Ad_g_minus * eta3

%% eta2_plus

exp_eta2_minus = [1 0 0 0;
                0 c2 -s2 0;
                0 s2 c2 0;
                0 0 0 1];

eta2 = [0 0 0 -1 0 0]';

tmp = gst_0_minus * exp_eta3_minus * exp_eta2_minus;

R = tmp(1:3, 1:3)

p = tmp(1:3, 4)

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R]

eta2_plus = Ad_g_minus * eta2

%% eta1_plus

exp_eta1_minus = [c1 s1 0 0;
                  -s1 c1 0 0;
                  0 0 1 0;
                  0 0 0 1];

eta1 = [0 0 0 0 0 1]';

tmp = gst_0_minus * exp_eta3_minus * exp_eta2_minus * exp_eta1_minus;

R = tmp(1:3, 1:3)

p = tmp(1:3, 4)

Ad_g_minus = [R, skew(p)*R;
              zeros(3, 3), R]

eta1_plus = Ad_g_minus * eta1
