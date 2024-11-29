%% Problème 1
clc
close all
clear all

load("APP6e_Formatif_Probleme1.mat")

Y1 = log(y1m);
X1 = xm;

C1 = [length(X1) sum(X1);
     sum(X1) sum(X1.^2)];
D1 = [sum(Y1);
      sum((Y1.*X1))];

A1 = pinv(C1)*D1;
b1 = A1(1);
m1 = A1(2);

Tau = -1/m1;
K = exp(b1);

Y2 = y2m;
X2 = xm.*sin(xm./7);

C2 = [length(X2) sum(X2);
     sum(X2) sum(X2.^2)];
D2 = [sum(Y2);
      sum((Y2.*X2))];

A2 = pinv(C2)*D2;
b2 = A2(1);
m2 = A2(2);

Alpha = m2;
Beta = b2;

Y1_appro = m1.*X1 + b1;
Y2_appro = m2.*X2 + b2;

y1_appro = K.*exp(-xm./Tau);
y2_appro = Alpha.*xm.*sin(xm./7) + Beta;

Y1_ = (1/length(Y1))*sum(Y1);
Y2_ = (1/length(Y2))*sum(Y2);

R2(1) = (sum((Y1_appro-Y1_).^2))/(sum((Y1-Y1_).^2));
R2(2) = (sum((Y2_appro-Y2_).^2))/(sum((Y2-Y2_).^2));

RMS(1) = sqrt((1/length(y1m))*sum((y1_appro - y1m).^2));
RMS(2) = sqrt((1/length(y2m))*sum((y2_appro - y2m).^2));

%% Probleme 2
clc 
close all
clear all

for n = 1:50
    x = n;
    iter = 0;
    tol = 10e-8;
        
    F1 = 80 * exp(-x/12);
    F2 = 3*x*sin(x/7) + 8;

    D1 = (-80/12)*exp(-x/12);
    D2 = (3*sin(x/7)) + ((3/7)*x*cos(x/7));

    F = F1 - F2;
    D = D1 - D2;

    while (abs(F) > tol && iter < 100)
        x = x - F/D;
       
        F1 = 80 * exp(-x/12);
        F2 = 3*x*sin(x/7) + 8;
    
        D1 = (-80/12)*exp(-x/12);
        D2 = (3*sin(x/7)) + ((3/7)*x*cos(x/7));
    
        F = F2 - F1;
        D = D2 - D1;

        iter = iter + 1;
    end
    
    Val(n, 1) = iter;
    Val(n, 2) = x;
    Val(n, 3) = F1;
    Val(n, 4) = F2;
end


%% Problème 4
clc
clear all
close all

tspan = [0 10];
z0 = [-10 -5 10];

options = odeset('AbsTol', 1e-6, 'RelTol', 1e-6);
[t, x] = ode45('probleme4', tspan, z0, options);

 figure(1)
 plot3(x(1,1), x(1,2), x(1,3),'go','Markersize',10) % début
 hold on
 plot3(x(:,1), x(:,2), x(:,3),'Linewidth',2) % trajectoire
 plot3(x(end,1), x(end,2), x(end,3),'rx','Markersize',10) % fin
 grid on
 xlabel('Position', 'Fontsize',15)
 ylabel('Vitesse' , 'Fontsize',15)
 zlabel('Accélération', 'Fontsize',15)
 legend('Début', 'Trajectoire', 'Fin')

%% Problème 5
clc
clear all
close all

x = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5];
F = [12.5 66.41 110 144.2 170 188.3 200 206.1 207.5 205.2 200 193.0 185 177.0 170 164.8 162.5 163.9 170 181.7 200];

figure 
grid on
hold on

h = x(2) - x(1);

N = length(F);

F1(1) = 0;
for n = 2:N
    F1(n) = (F(1) + F(n) + 2*sum(F(2:n-1))) * h/2;
end

plot(x, F1);

F2(1) = 0;
X2(1) = x(1);
for n = 3:2:N
    F2(((n-1)/2)+1) = (F(1) + F(n) + 4*sum(F(2:2:n-1)) + 2*sum(F(3:2:n-1)))*(h/3); 
    X2(((n-1)/2)+1) = x(n);
end

plot(X2, F2);