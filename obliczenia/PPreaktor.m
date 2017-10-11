LAMBDA = 10^-5;
S = 1938.9;
p_1 = -1.223e-04;
p_2 = -5.502e-05;
p_3 = -1.953e-04;
N0 = [10:10:100];

A = (p_1 .* N0) / LAMBDA;
B = (p_2 .* N0) / LAMBDA;
C = ((p_3 .* N0) / LAMBDA) + S;

delta = B.^2- 4.*A.*C;
Pdelta = sqrt(delta);

Q = (-B - Pdelta);
Q2 = (-B + Pdelta);
A1 = 2*A;

V_0 = Q./A1;
V_1 = Q2./A1;

