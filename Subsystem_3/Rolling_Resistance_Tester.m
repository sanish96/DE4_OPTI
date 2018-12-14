%% Wheels - Comparing Existing Product to Optimised Product

% Parameters
g = 9.81;
Wp = 100*g/2;
Wf = 50*g/2;
a = 0.015;

% Existing Product 
% w = 0.05;
% do = 0.22;
% di = 0.145;
% p = 1.25*10^3;
% E = 4.00*10^6;

% Optimised Product 
w = 0.04
do = 0.23
di = 0.2
p = 1.25*10^3;
E = 5.27*10^7;

Rollingresistance = (0.9.*(Wp+Wf+(p.*g*(((pi*((do/2)^2)*w)-(pi*((do-2*a)/2)^2)*(w-2*a))-(pi*((di/2)^2)*2*a)))).*sqrt((Wp+Wf+(p.*g*(((pi*((do/2)^2)*w)-(pi*((do-2*a)/2)^2)*(w-2*a))-(pi*((di/2)^2)*2*a))) )/E.*w))/do
