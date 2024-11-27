function f = code(t,z)
    u = 0;
    x_etoile = 10;
    
    Kp = 1;
    Kd = 1;
    Ka = 1;

    x0 = z(1);
    x1 = z(2);
    x2 = z(3);

    f1 = 1.25*x2 - 1.625*x1*abs(x1) - 1.75*x0^3;
    g1 = 1;
    ucmd = ((-(f1))+Kp*(x_etoile - x0) + Kd*(0-x1) + Ka*(0-x2))/g1;
    
    F = (f1) + g1*ucmd;

    f(1) = x1;
    f(2) = x2;
    f(3) = F;
    f(4) = ucmd > x_etoile;

    f = f(:);
end