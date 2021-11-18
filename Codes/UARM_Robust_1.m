clear all;
close all;
clc;

%%   Robot Parameters :

% Link Masses
m1 = 0.850;
m2 = 0.850;
m3 = 0.625;
% Link Inertias
I1 = 0.0075;
I2 = 0.0075;
I3 = 0.0060;
% Link Lengths
l1 = 0.203;
l2 = 0.203;
l3 = 0.203;
% Link COG Lengths
lc1 =0.096;
lc2 =0.096;
lc3 =0.077;

%% Computed Torque Controller Matrices (Phase 1):

Kp = [20];
Kv = [20];
Kp_n = rank(Kp);

%% Computed Torque Controller Matrices (Phase 2):

Kp = [10 0
      0  20];
Kv = [5  0
      0  20];
Kp_n = rank(Kp);
  
%% Performance Weighting Function Parameters (Phase 1):

Ms = 1.1;
omega_b = 1;
epsilon = 0.01;

%% Performance Weighting Function Parameters (Phase 2):

Ms = 1.5;
omega_b = 1;
epsilon = 0.01;

%% Control Weighting Function Parameters (Phase 1):

Mu = 20;
omega_bc = 10^6;
epsilon_1 = 1;

%% Control Weighting Function Parameters (Phase 2):

Mu = 50;
omega_bc = 10^6;
epsilon_1 = 1;

%% Weighting Functions State Space Realizations (Phase 1) :

W_e = tf([1/Ms omega_b],[1 omega_b*epsilon]);
W_delta = tf([1 omega_bc/Mu],[epsilon_1 omega_bc]);

sys_W_e = canon(W_e,'modal');
sys_W_delta = canon(W_delta,'modal');

A_e = sys_W_e.A;
B_e = sys_W_e.B;
C_e = sys_W_e.C;
D_e = sys_W_e.D;



A_delta = sys_W_delta.A;
B_delta = sys_W_delta.B;
C_delta = sys_W_delta.C;
D_delta = sys_W_delta.D;

%% Weighting Functions State Space Realizations (Phase 2) :
Numerator_e = {[1/Ms omega_b] [0];[0] [1/Ms omega_b]};
Denominator_e = [1 omega_b*epsilon];
W_e = tf(Numerator_e,Denominator_e);
Numerator_delta = {[1 omega_bc/Mu] [0];[0] [1 omega_bc/Mu]};
Denominator_delta = [epsilon_1 omega_bc];
W_delta = tf(Numerator_delta,Denominator_delta);

sys_W_e = canon(W_e,'modal');
sys_W_delta = canon(W_delta,'modal');

A_e = sys_W_e.A;
B_e = sys_W_e.B;
C_e = sys_W_e.C;
D_e = sys_W_e.D;

A_delta = sys_W_delta.A;
B_delta = sys_W_delta.B;
C_delta = sys_W_delta.C;
D_delta = sys_W_delta.D;

%% Augmented Plant State Space Matrices :

A = [zeros(Kp_n)            eye(Kp_n)           zeros(Kp_n)      zeros(Kp_n)
     -Kp                    -Kv                 zeros(Kp_n)      zeros(Kp_n)
     zeros(Kp_n)            zeros(Kp_n)         A_delta          zeros(Kp_n)
     B_e                    zeros(Kp_n)         zeros(Kp_n)      A_e];
 
B1 = [zeros(Kp_n)           zeros(Kp_n)
      eye(Kp_n)             zeros(Kp_n)
      zeros(Kp_n)           zeros(Kp_n)
      zeros(Kp_n)           B_e];
  
B2 = [zeros(Kp_n)
      eye(Kp_n)
      B_delta
      zeros(Kp_n)];
  
C1 = [zeros(Kp_n)           zeros(Kp_n)         C_delta         zeros(Kp_n)
      D_delta               zeros(Kp_n)         zeros(Kp_n)     C_e];
  
C2 = [eye(Kp_n)             zeros(Kp_n)         zeros(Kp_n)     zeros(Kp_n)];

D11 = [zeros(Kp_n)          zeros(Kp_n)
       zeros(Kp_n)          D_e];

D12 = [D_delta
       zeros(Kp_n)];
   
D21 = [zeros(Kp_n)          eye(Kp_n)];

D22 = [zeros(Kp_n)];

B = [B1 B2];
     
C = [C1;C2];

D = [D11 D12;D21 D22];

%Augmented Plant System
sys_aug = ss(A,B,C,D);

%Augmented Plant Matrix
P = sys_aug;
P_matrix = [A B;C D];
P_tf = tf(sys_aug);

%% H_infinity Controller Design :

[K_Hinf, sys_CL_Hinf, gamma, INFO_Hinf] = hinfsyn(P,Kp_n,Kp_n)
K_Hinf_tf = tf(K_Hinf)

figure(1);
bode(sys_CL_Hinf);

figure(2);
step(sys_CL_Hinf);

[U,S,V]=svd(sys_CL_Hinf.A);

figure(3);
loops = loopsens(ss(A,B2,C2,D22),K_Hinf); 
bode(loops.Si,'r',loops.Ti,'b',loops.Li,'g');
legend('S','T','Loop gain');

%% H2/Hinf Controller Design :

[K_Mix,sys_CL_Mix,normz_Hinf,INFO_Mix] = h2hinfsyn(P,Kp_n,Kp_n,Kp_n,[1,0])
K_Mix_tf = tf(K_Mix)

figure(1);
bode(sys_CL_Mix);

figure(2);
step(sys_CL_Mix);

figure(3);
loops = loopsens(ss(A,B2,C2,D22),K_Mix); 
bode(loops.Si,'r',loops.Ti,'b',loops.Li,'g');
legend('S','T','Loop gain');

%% H2 Controller Design :

[K_H2, sys_CL_H2, gamma_H2, INFO_H2] = h2syn(P,Kp_n,Kp_n)
K_H2_tf = tf(K_H2)

figure(1);
bode(sys_CL_H2);

figure(2);
step(sys_CL_H2);

figure(3);
loops = loopsens(ss(A,B2,C2,D22),K_H2); 
bode(loops.Si,'r',loops.Ti,'b',loops.Li,'g');
legend('S','T','Loop gain');

%% Uncertain System Model :
A_plant = [zeros(Kp_n) eye(Kp_n);-Kp -Kv];
B_plant = [zeros(Kp_n);eye(Kp_n)];
C_plant = [eye(Kp_n) zeros(Kp_n)];
D_plant = 0;

sys_plant = ss(A_plant,B_plant,C_plant,D_plant);
G0 = tf(sys_plant);
InputUnc = ultidyn('InputUnc',[Kp_n,Kp_n]);
Gpert = G0*(eye(Kp_n)+InputUnc*W_delta);
%P_uss = [W_e*Gpert;Gpert];
%P_uss = [W_e W_e*Gpert;eye(Kp_n) Gpert];
P_uss = [W_e;eye(Kp_n)]*[eye(Kp_n),Gpert];

%% Mu Controller Design :

figure(1);
fmu = logspace(-2,4,60);
opt = dksynOptions('FrequencyVector',fmu,'NumberofAutoIterations',5,'DisplayWhileAutoIter','on','MixedMU','on');
[K_DK, sys_CL_DK,bnd,INFO_DK] = dksyn(P_uss,Kp_n,Kp_n,opt)
K_DK_tf = tf(K_DK)

figure(2);
bode(sys_CL_DK);

figure(3);
step(sys_CL_DK);

figure(4);
loops = loopsens(ss(A,B2,C2,D22),K_DK); 
bode(loops.Si,'r',loops.Ti,'b',loops.Li,'g');
legend('S','T','Loop gain');
