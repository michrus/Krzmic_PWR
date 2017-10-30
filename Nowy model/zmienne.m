%Reaktor

B_1 = 0.000209;
B_2 = 0.001414;
B_3 = 0.001309;
B_4 = 0.002727;
B_5 = 0.000925;
B_6 = 0.000314;

B_t = 0.00689;

lambda_1 = 0.0125;
lambda_2 = 0.0308;
lambda_3 = 0.114;
lambda_4 = 0.307;
lambda_5 = 1.19;
lambda_6 = 3.19;

C_0 = 4625.71;
LAMBDA = 17.9 * 10^-6;
alfa_f = -1.1 * 10^-5;

alfa_c = -2.0 * 10^-4;
alfa_p = -1 * 10^-6;
M_ft = 222739;
C_PF = 0.059;
F_r = 0.974;

V_up = 1376;
V_lp = 1791;
V_hl = 1000;
V_cl = 2000;
V_mo = 540;
DENS_mo = 45.71;
C_PC = 1.39;

B_cl = -0.0631;
B_hl = -0.0741;

T_po0 = 539.5;
T_f10 = 1488.75;
T_mo10 = 548.303;
T_mo20 = 557.106;
T_f20 = 1506.38;
T_mo30 = 565.909;
T_mo40 = 574.711;
T_f30 = 1523.96;
T_mo50 = 583.514;
T_mo60 = 592.5;
T_hl0 = 592.5;
T_up0= 592.5;
T_lp0 = 539.5;
T_cl0 = 539.5;

ro_ex = 0;
ro_ex0 = 0;
POWER_i = 691244.4199;
T_set0 = 566;
tau_set = 30;
T_aves0 = 566;
tau_la1 = 10;
tau_la2 = 5;
tau_le = 80;
AUXCO_0 = 0;
T_clp0 = 539.5;
T_hlp0 = 592.5;
tau_rtd = 4;
React = 0.003375;
P_cor0 = 1;
P_cor10 = 3.436 * 10^9;
W_pth = 1.57559 * 10^8;

%obliczanie sta³ych parametrów dla reaktora

P_0 = (P_cor10 * 3.41) / (3600 * 3);
h = 200/3600;
A = 59900/3;
W_pt = W_pth/3600;

M_F = M_ft * 3;
M_C = V_mo * DENS_mo / 3;
M_up = V_up * DENS_mo;
M_lp = V_lp * DENS_mo;
M_hl = V_hl * DENS_mo;
M_cl = V_cl * DENS_mo;


tau_C = M_C / (2 * W_pt);
tau_up = M_up / W_pt;
tau_lp = M_lp / W_pt;
tau_hl = M_hl / W_pt;
tau_cl = M_cl / W_pt;


%dodatkowo do reaktrora na chwile

P = 849.7;
X8 = 1225.85;
K9 = -0.035;
H_st = P * K9 + X8;

Cl = 1.2195;
CP_fw = 1.218;
T_fi = 434.3;
POWER_s = P * Cl *(H_st - CP_fw * T_fi);

 
% Constant parameters of the UTSG MODEL
DENS_m = 530;
DENSW = 45.710;
DENS_r0 = 7.94;
DENS_d = 50.32;
DENS_dw = 47.66; 
DENS_g0 = 1.8325;
DENS_s = 52.32;
DENS_b0 = 13.614;

N = 3388;
D_o = 0.875;
D_i = 0.775;
L = 35.54;
L_s10 = 3.4517;
A_r = 48.7;
A_dw = 110.74;
A_d = 32;
A_fs = 60.67;
L_r = 9.63;
L_dw0 = 9.63;
L_d = 35.54;
V_p = 1077;
V_s = 3332.28;
V_r = 468.981;
V_dr = 4398.706;

T_h10 = 592.5;
T_pi0 = 592.5;
T_p10 = 587.36;
T_p20 = 557.4;
T_p40 = 539.5;
T_po0 = 539.5;
T_m10 = 553.49;
T_m20 = 536.129;
T_m40 = 527.38;
T_dw0 = 504.7;
T_d0 = 504.7;
T_sat0 = 521.9;
T_fi = 434.3;
T_fw = 434.3;
T_p30 = 542.5;
T_m30 = 529.38;

p_0 = 849.7;

H_f = 515.2;
H_fg = 678.3;
V_f = 0.02098;
V_fg = 0.5247;
X_e0 = 0.2

K1 = 3.5 * 10^-6;
K2 = -7.135 * 10^-4;
K3 = 0.17;
K4 = -0.2;
K5 = 0.14;
K6 = 0.14;
K7 = 2.37 * 10^-3;
K9 = -0.035;

X1 = 402.94;
X2 = 0.018;
X3 = 1.13096;
X4 = 370.751;
X5 = 850.04;
X6 = -0.181289;
X8 = 1225.85;

h_i = 1.25;
h_os = 0.87603;
h_ob = 1.87;
K_th = 0.0088275;

C_p1 = 1.39;
C_p2 = 1.165;
C_m = 0.11;
C_pfw = 1.2181;
W_fi0 = 1036.39;
W_pi = 10941.6;
W_10 = 5181.95;

C_l = 1.2195;
C_d = 4.101486234 * 10^-7;

TOU = 5;
TOU_1 = 250;
TOU_2 = 120;

G_1 = 65.2;
G_2 = 1.0;
G_v = 32.2;
W_nv = 0.63;
Z_tv = 3.18;

v0 = 0;
u0 = 0;
w0 = 0;
r0 = 0;
m0 = 0;

PI = 3.1415;

W_20 = W_10;
W_30 = W_10;
W_40 = W_10;
L_s20 = L - L_s10;
M_m = (DENS_m * N * L *PI * (D_o^2 - D_i^2)) / (4*144);
M_m1 = (M_m * L_s10) / L;
M_m4 = M_m1;
M_m2 = M_m * (L_s20 / L);
M_m3 = M_m2;
S_m = (L * PI * D_o * N)/12;
S_ms1 = S_m * (L_s20 / L);
S_ms3 = S_ms2;
S_ms4 = S_ms1;
S_pm1 = (S_ms1 * D_i) / D_o;
S_pm2 = (S_ms2 * D_i) / D_o;
S_pm3 = S_pm1;
S_pm4 = S_pm2;
P_r1 = S_pm1 / L_s10;
P_r2 = S_pm2 / L_s20;
D_m = (D_i + D_o) / 2;

A_p = (PI * (D_i^2) * N) / (4 * 144);
M_p = (DENS_w * A_p * L);
M_p1 = (M_p * L_s10) / L;
M_p2 = (M_p * L_s20) / L;
M_p3 = M_p2;
M_p4 = M_p1;
M_s1 = A_fs * DENS_s * L_s10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_pm = 1 / (((1 / h_i) + (D_i * log(D_m/D_i))) / (24 * K_th));
U_ms1 = 1 / (((1 / h_os) + (D_o * log(D_o/D_m))) / (24 * K_th));
U_ms2 = 1 / (((1 / h_ob) + (D_o * log(D_o/D_m))) / (24 * K_th));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_pi = (V_p - A_p * 2 * L) / 2;
M_pi = DENS_w * V_pi
M_po = M_pi;
M_d = DENS_d * A_d *L_d;

tau_pi = M_pi / W_pi;

H_b0 = H_f + ((X_e0 * H_fg) / 2);
H_xe0 = H_f + (X_e0 * H_fg);
DENS_b0 = 1/(V_f + (X_e0 * V_fg) / 2) 
L_b0 = L_s20;

C_1 = 1 / sqrt(C_d);
W_r = (1 - X_e0) * W_40;

DEN_1 = (A_fs * DENS_s * C_p2 * (T_d + X1 + (K5 * P))/2);
DEN_2 = (DENS_b * A_fs * (L - L_s1) * ((X5 + K4 * P) / 2));

K01 = -(K1 + ((K2 * X_e) / 2



















