%% Optimisation

clc

% Paramaters 
g = 9.81;
Wp = 100*g/2;
Wf = 150*g/2;

% Variables (Set to constants equal to current specifications)
w = 0.05;
do = 0.22;
di = 0.145;

% Test dimensional Values

p_upper = 1.8*10^3;
p_lower = 0.92*10^3;

E_lower = 2.67*10^6;
E_upper = 5.27*10^7;

%Averages
p = 1.07071*10^3;
E = 1.39*10^7;

%% Plotting 

subplot(1,2,1)

p_values = linspace(p_lower, p_upper, 100);
R = (0.9.*(Wp+Wf+p_values.*pi.*w.*g*((do^2-di^2)./4)).*sqrt((Wp+Wf+p_values.*pi.*w.*g.*((do^2-di^2)./4))./E.*w))/do
plot(p_values,R)
title('Effect of Density on Rolling Resistance')
xlabel('Density (p)')
ylabel('Rolling Resistance (N)')
set(gca,'Color','w')

subplot(1,2,2)

E_values = linspace(E_lower, E_upper, 100);
R2 = (0.9.*(Wp+Wf+p.*pi.*w.*g*((do^2-di^2)./4)).*sqrt((Wp+Wf+p.*pi.*w.*g.*((do^2-di^2)./4))./E_values.*w))/do
plot(E_values,R2)
title('Effect of Elastic Modulus on Rolling Resistance')
xlabel('Elastic Modulus (E)')
ylabel('Rolling Resistance (N)')
set(gca,'Color','w')

