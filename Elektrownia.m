% dN = ((P_1*V^2 + P_2*V + P_3) / LAMBDA)* N + S;
% 
% dMpc = M_in - M_out;
% 
% dTpc = (1/(Cp_pc * M_pc)) * (Cp_pc * M_in*(T_pc_i - T_pc) + C_p_pc * M_out(T_pc - T_pc + 15) + W_r - W_sg_6 - W_loss_pc);
% 
% dTsg = (1 / C_lp_sg * M_sg) * (C_lp_sg * m_sg * T_sg_sw - C_vp_sg * T_sg - m_sg * E_evap_sg + Kt_sg * (T_pc - T_sg) - W_loss_sg);
% 
% dTpr = (1 / C_p_pr * M_pr) * (X_mpr_gtz * C_p_pc * M_pr * T_pc_hl + X_mpr_ltz * C_p_pr * M_pr * T_pr - Cp_pr * M_pr * T_pr - W_loss_pr + W_heat_pr);
% 
% P_sg = 28884.78 - 258.01 * T_sg + 0.63 * T_sg^2;
% 
% PSI_N = C_psi * N;
% 
% JOTA_pr = (1 \ A_pr) * ((M_pc/(FI_pc_Tpc));
% 
% FI_pc_Tpc = C_fi_0 + C_fi_1 * T_pc + C_fi_2 * T_pc ^ 2;


C_fi_0 = 581.2;
C_fi_1 = 2.98;
C_fi_2 = -0.008448;



% tabela 2 [1]
P_1 = -1.223 * (10 ^ -4);
P_2 = -5.502 * (10 ^ -5);
P_3 = -1.953 * (10 ^ -4);
S = 1938.9;

% tabela 3 [1]

C_p_pc = 5355;
K_t_sg = 9.5296 * (10 ^ 6);
W_loss_pc =  2.996 * (10 ^ 7);

% tabela 4 [1]

M_sg = 34920;
W_loss_sg = 1.8932  * (10 ^ 7);
C_vp_sg = 3635.6;
C_lp_sg = 3809.9;

% tabela 5

C_p_pr = 6873.1;
W_loss_pr = 1.6823 * (10 ^ 5);

%tabela 6 [1]
N = 99.3;
V = 0;
M_pc = 200000;
T_pc = 281.13;
T_pc_hl = 296.13;
T_pc_cl = 266.13;
T_pc_i = 258.85;
M_out = 2.11;
M_in = 1.4222;
% JOTA_pr = 4.8;
M_pr = 19400;
T_pr = 326.57;
P_pr = 123;
W_heat = 168;
M_sg = 34920;
% JOTA_sg = 1.85;
T_sg = 257.78;
P_sg = 45.3;
M_sg_ss = 119.31;
M_sg_sw = 119.31;
m_sg = 119.31;
T_sg_sw = 220.85;
W_r = 13.75 * (10 ^ 8);
W_sg_6 = 13.351 * (10 ^ 8);


% tabela 7 [1]

C_psi = 13.75 * (10 ^ 6);
B = 0.0064;
lambda = 0.1;
LAMBDA = 10 ^ -5;
V_0_pc = 242;
M_0_pc = 180600;
A_pr = 4.52;
V_pr_vessel = 44;
E_evap_sg = 1.658 * (10 ^ 6);
Kt_sg = 9.5296*(10^6);

P_1 = -1.223e-4;
P_2 = -5.502e-5;
P_3 = -1.953e-4;
S = 1938.9;

