% plot(Time,N_out);
% grid on;
% axis([10 100 -10 100]);
% 




K = 1;
T = 15.79;
T0 = 3.52;

%0%
%   P
Kp1 = (0.3*T)/(K*T0);

%   PI
Kp2 = (0.35*T)/(K*T0);
Ti2 = (1.2*T0);

%   PID


Kp3 = (0.6*T)/(K*T0);
Ti3 = T;
Td3 = (0.5*T0); 