%% Probleme 1
clc 
close all
clear all

load("APP6e_Formatif_Probleme1.mat");

Y1 = log(y1m);
X1 = xm;

C1 = [length(X1) sum(X1);
      sum(X1) sum(X1.^2)];
D1 = [sum(Y1);
      sum(Y1.*X1)];

A1 = pinv(C1)*D1;

Y2 = y2m;
X2 = xm.*sin(xm/7);

C2 = [length(X2) sum(X2);
      sum(X2) sum(X2.^2)];
D2 = [sum(Y2);
      sum(Y2.*X2)];

A2 = pinv(C2)*D2;

Tau = -1/A1(2);
K = exp(A1(1));
Alpha = A2(2);
Beta = A2(1);


y1_appro = K*exp(-xm/Tau);
y2_appro = Alpha.*xm.*sin(xm/7) + Beta;
Y1_appro = (-1/Tau).*X1 + log(K);
Y2_appro = Alpha*X2 + Beta;

Y1_ = (1/length(Y1))*sum(Y1);
Y2_ = (1/length(Y2))*sum(Y2);

R2(1) = sum((Y1_appro-Y1_).^2) / sum((Y1-Y1_).^2);
R2(2) = sum((Y2_appro-Y2_).^2) / sum((Y2-Y2_).^2);

RMS_abs(1) = sqrt((1/length(y1m))*sum((y1_appro-y1m).^2));
RMS_abs(2) = sqrt((1/length(y2m))*sum((y2_appro-y2m).^2));

RMS_rel(1) = sqrt((1/length(y1m))*sum(((y1_appro-y1m).^2)./y1m));
RMS_rel(2) = sqrt((1/length(y2m))*sum(((y2_appro-y2m).^2)./y2m));

%% Problème 2
clc
close all
clear all

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


%% Problème 4
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
Erreur = ((h^2)/12)*(Fpb-Fpa)

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


% b)
m = 1687.5;
Vitesse = sqrt(F_appro1(end)/m);
Vitesse = Vitesse*(3600/1000)

%% Problème 6
clc
clear all
close all

% a)
% 3^5

% b)
% Tu l'essaye en regardant avec un gros et un petit si les graphiques ce
% ressemblent cela veux dire que tu as une erreur très petite

% c)
% Trop petit prends trop de temps de calcul
% Trop grand pas assez de précision

%% Problème 7
clc
clear all
close all

load("APP6e_Formatif_Probleme8.mat")

h = t(2) - t(1);

% a)
% Trapèze
F(1) = 0; % Initialisation de l’intégrale
for n = 2:length(acc_mes)
    F(n) =  (acc_mes(1) + acc_mes(n) + 2*sum(acc_mes(2:n-1))) * h/2;
end 

%Erreur Trapèze
Fpa = (acc_mes(2)-acc_mes(1))/h;
Fpb = (acc_mes(end)-acc_mes(end-1))/h;
Erreur = ((h^2)/12)*(Fpb-Fpa)

% b)
%Simpson
Fs(1) = 0;
Xs(1) = t(1);
for n = 3:2:length(acc_mes)
    Fs(((n-1)/2)+1) = (acc_mes(1) + acc_mes(n) + (4*sum(acc_mes(2:2:n-1))) + (2*sum(acc_mes(3:2:n-1))))*(h/3);
    Xs(((n-1)/2)+1) = t(n);
end

%Erreur Simpsons
Fpppa = (acc_mes(4) - 3*acc_mes(3) + 3*acc_mes(2) - acc_mes(1) )/h^3;
Fpppb = (acc_mes(end) - 3*acc_mes(end-1) + 3*acc_mes(end-2) - acc_mes(end-3) )/h^3;
Erreur = ((h^4)/180)*(Fpppb-Fpppa)

% figure
% hold on
% plot(Xs, Fs);
% plot(t, F)


%% Problème 10
clc
clear all
close all

tspan = [0 10];
z0 = [2 -0.1 -0.1];

options = odeset('AbsTol', 1e-6, 'RelTol',1e-6);
[t, x] = ode45('Probleme10', tspan, z0, options);

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