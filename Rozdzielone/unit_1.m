c_phi_0 = 581.2;
c_phi_1 = 2.98;
c_phi_2 = -0.008448;



% tabela 2 [1]
P_1 = -1.322 * (10 ^ -4);
P_2 = -6.08 * (10 ^ -5);
P_3 = -2.85 * (10 ^ -4);
S = 2859;

% tabela 3 [1]

c_p_PC = 5415;
K_T_SG = 8.4004 * (10 ^ 6);
W_loss_PC =  2.2469 * (10 ^ 8);

% tabela 4 [1]

M_SG = 35611;
W_loss_SG = 1.9166  * (10 ^ 5);
cv_p_SG = 3489;
cl_p_SG = 3871.3;

% tabela 5

c_p_PR = 6873.1;
W_loss_pr = 1.6823 * (10 ^ 5);

%tabela 6 [1]
N = 100.3;
V = -0.0025;
M_pc = 200000;
T_PC = 278.03;
T_PC_HL = 293.08;
T_pc_cl = 262.98;
T_PC_I = 246.1;
m_out = 2.9722;
% m_out = 1.4222;
m_in = 1.4222;
% JOTA_pr = 4.8;
M_pr = 19400;
T_pr = 326.57;
P_pr = 123;
W_heat = 168;
M_SG = 35611;
% JOTA_sg = 1.85;
T_SG = 255.13;
P_sg = 43.3;
M_sg_ss = 120.56;
M_sg_sw = 120.56;
m_SG = 120.56;
m_sg = 119.31;
T_SG_SW = 219.65;
W_r = 13.75 * (10 ^ 8);
W_SG = 11.542 * (10 ^ 8);


% tabela 7 [1]

C_psi = 13.75 * (10 ^ 6);
B = 0.0064;
lambda = 0.1;
LAMBDA = 10 ^ -5;
V0_PC = 242;
M_0_pc = 180600;
A_pr = 4.52;
V_pr_vessel = 44;
E_evap_SG = 1.658 * (10 ^ 6);


% zmienne wej�ciowe PR

W_heat_PR = 1; 
W_loss_PR = 1.6823 * 10^5;

M_PR = M_pc-V0_PC * (581.2 + 2.98*T_PC - 0.0848*(T_PC^2));  
A_PR = 4.52
