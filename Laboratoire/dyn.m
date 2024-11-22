function f = dyn(t,z)

    x1 = z(1);
    x2 = z(2);

    ut = 0.0;

    f(1) = x1;
    f(2) = ut -0.6*x2 - 3*x1 - x1^2;
    
    f = f(:);

end

