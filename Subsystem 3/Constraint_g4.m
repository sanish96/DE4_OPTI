%% Constraint g4: Inner diameter is less than outer diameter

function g3 = Constraint_g3(x) 
    % Parameters 
    t = 0.015;

    % Variables
    w = x(1);
    
    % Constraint g1 (dimensional)
    g3 = 0;
    