%% Wheels - Comparing Mass To Existing Product

% parameters
p = 1250;
a = 0.015;

% existing dimensions
% w = 0.05;
% do = 0.22;
% di = 0.145;

% my dimensions
w = 0.04
do = 0.23
di = 0.2

% function
v = (((pi*((do/2)^2)*w)-(pi*((do-2*a)/2)^2)*(w-2*a))-(pi*((di/2)^2)*2*a));
mass = v*p