function [c, ceq] = nlcon(x)
k = 15.5;
rho = 7950;
E = 197e6;
tb = 5;
T_amb = 25;

Fb = 2994.8;
Q = 8984.4;
Ap = pi*0.01^2;
rho = 7950;

% Equations
A = (pi * x(1)^2) - ((pi * x(4))/72) * (x(2)^2 - x(3)^2);
T_max = ((Q * x(5))/(tb * k * A)) + T_amb;
delta = (Fb * x(1))/(A * E);
sigma_max = (6 * Fb * x(1))/((x(1) * sin(x(4)/2)) * x(5)^2);

% Inequalities
c(1) = T_max - 250;
c(2) = delta - 0.005;
c(3) = sigma_max - 50e6;

ceq = [];