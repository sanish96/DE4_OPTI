%% Constraint g2: Mass is greater than zero

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
    