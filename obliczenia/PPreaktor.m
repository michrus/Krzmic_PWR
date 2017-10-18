LAMBDA = 10^-5;
S = 2859;
p_1 = -1.322 * (10 ^ -4);
p_2 = -6.08 * (10 ^ -5);
p_3 = -2.85 * (10 ^ -4);
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

