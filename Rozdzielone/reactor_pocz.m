C_psi = 13.75 * 10^6;
P_1 = -1.223e-04;
P_2 = -5.502e-5;
P_3 = -1.953e-04;
V = 0;

S = 1938.9;
lambda = 0.1;
LAMBDA = 10^-5;
B = 0.0064;

N_roz = -(S)./((((P_1.*V.^2 + P_2.*V +P_3)-B)./LAMBDA) + ((B)./(LAMBDA)));

N_pod = (-S)./((P_1.*V.^2 + P_2.*V +P_3)./LAMBDA);

%gowno
%eksperymenty eksperymenty