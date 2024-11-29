function f = probleme4(t,z) 
    x_etoile = 10;
    
    x0 = z(1);
    x1 = z(2);
    x2 = z(3);
    
    f1 = -1.75*x2 - 1.625*x1*abs(x1) - 1.25*x0^3;
    g1 = 1;

    Kd = 28;
    Ka = 9;
    Kp = 40;
    
    
    ucmd = (-f1 + Kp*(x_etoile - x0) + Kd*(0-x1) + Ka*(0-x2))/g1;
    
    f(1) = x1;
    f(2) = x2;
    f(3) = f1 + g1*ucmd;

    f = f(:);
end