function f = apollo(t,z)
  u = 1.0/(82.45);

  x = z(1);
  y = z(2);

  vx = z(3);
  vy = z(4);

  u_etoile = 1 - u;

  r1 = sqrt(((x+u)^2)+y^2);
  r2 = sqrt(((x-u_etoile)^2)+y^2);

  f(1) = vx;
  f(2) = vy;
  f(3) = (2*vy) + x - ((u_etoile*(x+u))/(r1^3)) -((u*(x-u_etoile))/(r2^3));
  f(4) = (-2*vx) + y - ((u_etoile*y)/(r1^3)) -((u*y)/(r2^3)) ;
  
  f = f(:);
