function [c,ceq] = nlcon(x)
    c = [Constraint_g1(x); Constraint_g2(x); Constraint_g3(x)];
    ceq = 0;
end 