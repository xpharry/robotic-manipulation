%% howework 4

%% Q3

% variables
syms c1 s1 c2 s2 c3 s3 c4 s4 c5 s5 c6 s6
syms l0 l1 l2

% axis
w1 = [0, 0, 1]';
w2 = [0, 1, 0]';
w3 = [-1, 0, 0]';
w4 = [-1, 0, 0]';
w5 = [-1, 0, 0]';
w6 = [0, 1, 0]';

% rotation matrix
Rz1 = [c1, -s1, 0;
       s1, c1, 0;
       0, 0, 1];   
Ry2 = [c2, 0, -s2;
       0, 1, 0;
       s2, 0, c2];  
R_x3 = [1, 0, 0;
        0, c3, s3;
        0, -s3, c3]; 
R_x4 = [1, 0, 0;
        0, c4, s4;
        0, -s4, c4];
R_x5 = [1, 0, 0;
        0, c5, s5;
        0, -s5, c5];

% new axis after rotation
w2_ = Rz1 * w2;
w3_ = Rz1 * Ry2 * w3;
w4_ = Rz1 * Ry2 * R_x3 * w4;
w5_ = Rz1 * Ry2 * R_x3 * R_x4 * w5;
w6_ = Rz1 * Ry2 * R_x3 * R_x4 * R_x5 * w6;

% new origins after rotation
q1 = [0; 0; l0];
q2 = [0; 0; l0];
q3 = [0; 0; l0];
q4_ = q1 + Rz1 * Ry2 * R_x3 * [0; l1; 0];
q5_ = q4_ + R_x4 * [0; l2; 0];
q6_ = q5_;

% construct eta
eta1 = [-cross(w1, q1); w1];
eta2 = [-cross(w2_, q2); w2_];
eta3 = [-cross(w3_, q3); w3_];
eta4 = [-cross(w4_, q4_); w4_];
eta5 = [-cross(w5_, q5_); w5_];
eta6 = [-cross(w6_, q6_); w6_];

% jacobian
Jst = [eta1, eta2, eta3, eta4, eta5, eta6];
disp('Jst = ')
disp(Jst)

