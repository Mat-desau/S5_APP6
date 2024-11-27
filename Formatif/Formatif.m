clc
close all
clear all

%% Problème 2
for n = 1:50
    x = n;

    tol = 10e-8;
    iter = 0;

    F1 = 80*exp(-x/12);
    D1 = (-80/12)*exp(-x/12);
    
    F2 = 3*x*sin(x/7) + 8;
    D2 = 3*sin(x/7) + ((3*x)/7)*cos(x/7);
    
    F = F1 - F2;
    D = D1 - D2;
    
    while (abs(F) > tol && iter < 100)
        x = x - F/D;
    
        F1 = 80*exp(-x/12);
        D1 = (-80/12)*exp(-x/12);
        
        F2 = 3*x*sin(x/7) + 8;
        D2 = 3*sin(x/7) + ((3*x)/7)*cos(x/7);
        
        F = F1 - F2;
        D = D1 - D2;
    
        iter = iter + 1;
    end
    Val(n,1) = iter;
    Val(n,2) = x;
end


%% Problème 3
clc 
close all
clear all

tspan = [0 10];
z0 = [-10 -5 10 0];

options = odeset('abstol', 1e-6, 'RelTol', 1e-6);
[t, x] = ode45('code', tspan, z0, options);

% figure(1)
% plot3(x(1,1), x(1,2), x(1,3),'go','Markersize',10) % début
% hold on
% plot3(x(:,1), x(:,2), x(:,3),'Linewidth',2) % trajectoire
% plot3(x(end,1), x(end,2), x(end,3),'rx','Markersize',10) % fin
% grid on
% xlabel('Position', 'Fontsize',15)
% ylabel('Vitesse' , 'Fontsize',15)
% zlabel('Accélération', 'Fontsize',15)
% legend('Début', 'Trajectoire', 'Fin')

figure;
plot(t, x(:,1));

%% Problème 5
clc
clear all
close all

x = [0 0.025 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.225 0.25 0.275 0.3 0.325 0.35 0.375 0.4 0.425 0.45 0.475 0.5];
F = [12.5 66.41 110 144.2 170 188.3 200 206.1 207.5 205.2 200 193.0 185 177.0 170 164.8 162.5 163.9 170 181.7 200];

figure;
hold on

h = x(2) - x(1);

% Trapèze
F_appro1(1) = 0; % Initialisation de l’intégrale
for n = 2:length(F)
    F_appro1(n) =  (F(1) + F(n) + 2*sum(F(2:n-1))) * h/2;
end 

%Erreur Trapèze
Fpa = (F(2)-F(1))/h;
Fpb = (F(end)-F(end-1))/h;
Erreur = (h^2)/(Fpb-Fpa)

plot(x, F_appro1)

%Simpson
Fs(1) = 0;
Xs(1) = x(1);
for n = 3:2:length(F)
    Fs(((n-1)/2)+1) = (F(1) + F(n) + (4*sum(F(2:2:n-1))) + (2*sum(F(3:2:n-1))))*(h/3);
    Xs(((n-1)/2)+1) = x(n);
end

plot(Xs, Fs)

%Erreur Simpsons
Fpppa = (F(4) - 3*F(3) + 3*F(2) - F(1) )/h^3;
Fpppb = (F(end) - 3*F(end-1) + 3*F(end-2) - F(end-3) )/h^3;
Erreur = ((h^4)/180)*(Fpppb-Fpppa)

%% Problème 6
clc
clear all
close all

