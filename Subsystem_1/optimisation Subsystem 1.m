

%% Parameters

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

%% Optimisation

objective = @(x) ((pi * x(1)^2) - ((pi * x(4))/72) * (x(2)^2 - x(3)^2)) * x(5) * rho;

% variable bounds
lb = [0.04, 0.02, 0.02, 0, 0.002];
ub = [0.11, 0.09, 0.09, 72, 0.004];

% initial guess
% rng(0,'twister');
% r = [(ub(1)-lb(1)).*rand(1000,1) + lb(1), (ub(2)-lb(2)).*rand(1000,1) + lb(2), (ub(3)-lb(3)).*rand(1000,1) + lb(3), (ub(4)-lb(4)).*rand(1000,1) + lb(4), (ub(5)-lb(5)).*rand(1000,1) + lb(5)];
% y = linspace(lb(1),ub(1),100);

% for row = 1:1000
%     x0 = r(row, :);
%     
%     A = [0,-1,1,0,0; -1,1,0,0,0];
%     b = [0; -0.02];
%     Aeq = [];
%     beq = [];
%     
%     nonlincon = @nlcon;
%     
%     x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon);
    
    
    
% x0 = [0.11, 0.09, 0.09, 72, 0.004];
x0 = [0.8, 0.04, 0.02, 10, 0.002];
% % show initial objective
% disp(['Initial Objective: ' num2str(objective(x0))])

% linear constraints
A = [0,-1,1,0,0; -1,1,0,0,0];
b = [0; -0.02];
Aeq = [];
beq = [];

% % nonlinear constraints
nonlincon = @nlcon;

% optimize with fmincon
%[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] 
% = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon);

% rng default % For reproducibility
% gs = GlobalSearch;
% sixmin = objective;
% problem = createOptimProblem('fmincon','x0',x0, 'objective',sixmin,'lb',lb,'ub',ub,'nonlcon', nonlincon);
% x = run(gs,problem);

% show final objective
disp(['Final Objective: ' num2str(objective(x))])

% print solution
disp('Solution')
disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])
disp(['x4 = ' num2str(x(4))])
disp(['x5 = ' num2str(x(5))])