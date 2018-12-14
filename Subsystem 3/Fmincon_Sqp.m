
%% Optimisation Using Fmincon Sqp.

% Measuring run time
tic

% Parameters 
g = 9.81;
Wp = 100*g/2;
Wf = 50*g/2;
p = 1.25*10^3;
E = 5.27*10^7
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
nonlincon = @nlcon;

% Objective function
objective = @(x) (0.9.*(Wp+Wf+(p.*g*(((pi*((x(2)/2)^2)*x(1))-(pi*((x(2)-2*t)/2)^2)*(x(1)-2*t))-(pi*((x(3)/2)^2)*2*t)))).*sqrt((Wp+Wf+(p.*g*(((pi*((x(2)/2)^2)*x(1))-(pi*((x(2)-2*t)/2)^2)*(x(1)-2*t))-(pi*((x(3)/2)^2)*2*t))) )/E.*x(1)))/x(2)

% Optimisation with fmincon sqp
options = optimoptions('fmincon', 'Display', 'iter' , 'Algorithm','sqp');
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon, options);

% Display solution
disp('Solution')
disp(['w = ' num2str(x(1))])
disp(['do = ' num2str(x(2))])
disp(['di = ' num2str(x(3))])

toc
%% Functions 

function [c,ceq] = nlcon(x)
    c = [Constraint_g1(x); Constraint_g2(x)];
    ceq = 0;
end 








