function f = Probleme10(t,z)
    u = 0;

    x0 = z(1);
    x1 = z(2);
    x2 = z(3);
    
    x_etoile = 4;

    Kp = 4;
    Kd = 8;
    Ka = 3;
    
    f1 = -2*x2 - 2*x1 - 20*x1*abs(x1) - 4*x0 + x0^3;
    g1 = 1;
    ucmd = ((-(f1))+Kp*(x_etoile - x0) + Kd*(0-x1) + Ka*(0-x2))/g1;
    
    F = (f1) + g1*ucmd;

    f(1) = x1;
    f(2) = x2;
    f(3) = F;
    % f(3) = ucmd - (2*x2) - (2*x1) - (20*x1*abs(x1)) - (4*x0) + (x0^3);

    f = f(:);
end