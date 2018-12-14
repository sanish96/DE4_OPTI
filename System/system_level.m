clear variables

%% System Parameters

battery_mass = 2;
handlebar_mass = 0.5;
steering_column_mass = 2;
front_fork_mass = 1.5; 

%% Deck Optimisation Subsystem 2 for minimum mass result

deck_mass = 3.7; % mass of deck result from Subsystem 2 Optimisation

%% Subsytem 1 Wheel Reoptimisation Using Fmincon Interior Points

% Measuring run time 
tic 

% Parameters 
g = 9.81;
Wp = 100*g/2;
Wf = deck_mass*g/2; % mass of deck result from Subsystem 2 Optimisation
p = 1.25*10^3;
E = 5.27*10^7;
t = 0.015;

% Initial guess
x0 = [0.05,0.22,0.145];

% Variable bounds
lb = [0.04 0.1 0.08];
ub = [0.08 0.25 0.20];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Nonlinear constraints
nonlincon = @nlcon3;

% Objective function
objective3 = @(x) (0.9.*(Wp+Wf+(p.*g*(((pi*((x(2)/2)^2)*x(1))-(pi*((x(2)-2*t)/2)^2)*(x(1)-2*t))-(pi*((x(3)/2)^2)*2*t)))).*sqrt((Wp+Wf+(p.*g*(((pi*((x(2)/2)^2)*x(1))-(pi*((x(2)-2*t)/2)^2)*(x(1)-2*t))-(pi*((x(3)/2)^2)*2*t))) )/E.*x(1)))/x(2);

% Optimisation with fmincon interior points
x = fmincon(objective3,x0,A,b,Aeq,beq,lb,ub,nonlincon);

% % Display solution
% disp('Solution')
% disp(['w = ' num2str(x(1))])
% disp(['do = ' num2str(x(2))])
% disp(['di = ' num2str(x(3))])

toc

%% Brake Reoptimisation Using Fmincon Interior Points


% Known parameters
k = 15.5;
rho = 7950;
E = 197e6;
tb = 5;
T_amb = 25;

% Calculated parameters
Fb = 2994.8;
Q = 8984.4;
Ap = pi*0.01^2;
rho = 7950;


%Optimisation

objective1 = @(z) ((pi * z(1)^2) - ((pi * z(4))/72) * (z(2)^2 - z(3)^2)) * z(5) * rho;

% variable bounds
%VARIABLE 1 BOUNDS ARE FIXED TO VALUE FROM SUBSYSTEM 2
lb = [x(3)/2, 0.02, 0.02, 0, 0.002]; 
ub = [x(3)/2, 0.09, 0.09, 72, 0.004];

% initial guess
z0 = [0.09, 0.07, 0.07, 62, 0.002];

% linear constraints
A = [0,-1,1,0,0; -1,1,0,0,0];
b = [0; -0.02];
Aeq = [];
beq = [];

% % nonlinear constraints
nonlincon = @nlcon1;

% optimize with fmincon
options = optimoptions(@fmincon,'Algorithm','interior-point');
[z, brake_mass] = fmincon(objective1,z0,A,b,Aeq,beq,lb,ub,nonlincon,options);

% % show final objective
% disp(['Final Objective: ' num2str(objective1(z))])
% 
% % print solution
% disp('Solution')
% disp(['x1 = ' num2str(z(1))])
% disp(['x2 = ' num2str(z(2))])
% disp(['x3 = ' num2str(z(3))])
% disp(['x4 = ' num2str(z(4))])

%% Mass of Brake

mass_brake = brake_mass;

%% Mass of Deck

mass_deck = deck_mass;

%% Mass of Wheel

% parameters
p = 1250;
a = 0.015;

% my dimensions
w = x(1);
do = floor(x(2)*100)/100;
di = x(3);

% function
v = (((pi*((do/2)^2)*w)-(pi*((do-2*a)/2)^2)*(w-2*a))-(pi*((di/2)^2)*2*a));
mass_wheel = v*p;

%% Mass of System

system_mass = mass_brake + mass_deck + 2 * mass_wheel + battery_mass + handlebar_mass + steering_column_mass + front_fork_mass;

disp(['System Mass(kg) = ' num2str(system_mass)])

%% Functions 

% Subsystem 1 Brake

function [c, ceq] = nlcon1(z)
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
A = (pi * z(1)^2) - ((pi * z(4))/72) * (z(2)^2 - z(3)^2);
T_max = ((Q * z(5))/(tb * k * A)) + T_amb;
delta = (Fb * z(1))/(A * E);
sigma_max = (6 * Fb * z(1))/((z(1) * sin(z(4)/2)) * z(5)^2);

% Inequalities
c(1) = T_max - 250;
c(2) = delta - 0.005;
c(3) = sigma_max - 50e6;

ceq = [];

end

% Subsystem 3 Wheel

function [c,ceq] = nlcon3(x)
    c = [Constraint_g1(x); Constraint_g2(x)];
    ceq = 0;
end 

% Subsystem 3 Wheel - Constraint g2: Mass is greater than zero

function g2 = Constraint_g2(x) 
    % Parameters 
    p = 1.25*10^3;
    t = 0.015;
    
    %Variables
    w = x(1);
    do = x(2);
    di = x(3);
    
    % Constraint g1 (mass)
    g2 = -p*(((pi*((do/2)^2)*w)-(pi*((do-2*t)/2)^2)*(w-2*t))-(pi*((di/2)^2)*2*t)) + 0;
end
    
% Subsystem 3 Wheel - Constraint g1: Mass is greater than current model

function g1 = Constraint_g1(x)
    % Parameters 
    p = 1.25*10^3;
    t = 0.015;
    
    % Variables
    w = x(1);
    do = x(2);
    di = x(3);
    
    % Constraint g1 (mass)
    g1 = p*(((pi*((do/2)^2)*w)-(pi*((do-2*t)/2)^2)*(w-2*t))-(pi*((di/2)^2)*2*t)) - 0.6;
end
    
