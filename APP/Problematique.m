clc
clear all
close all
warning off


%% Assignation des variables
m = 50; %kg
J = 1.5; %kg-m^2
R_mars = 3397e03; %m
U_mars = 42830e09; %m^3/^2
S = 0.8; %s^2
d = 0.05; %m
C_D0 = 1.2;
C_Lalpha = 0.80;
C_Malpha = -0.07;
C_Mq = -0.05;
C_Mgamma = 0.10;

%% Conditions initiales
V_ini = 6100; %m/s
Gamma_ini = -20.5; %degré
h_ini = 120000; %m
s_ini = 0.0; %degré
Theta_ini = -80; %deg
q_ini = 0.0; %degré/s

%% Conditions finales désirer
V_fin1 = 250; %m/s
V_fin2 = 300; %m/s
h_fin = 10000; %m

%% Conditions initiales NASA
V_ini_NASA = 6100; %m/s
Gamma_ini_NASA = -90; %degrée
H_ini_NASA = 120000; %m
s_ini_NASA = 0.0; %degrée
Theta_ini_NASA = -90; %degrée
q_ini_NASA = 0.0; %degrée

%% Load les valeurs
load("Accelero_Data_from_NASA.mat")

%% Accélération
Accelertion = struct("X", t, "Y", acc_mes);

%% Méthode de trapèze
Vitesse_mes(1,1) = V_ini_NASA;
h = t(2);
for n = 2:length(acc_mes)
    Vitesse_mes(n,1) = Vitesse_mes(1,1) - (acc_mes(1) + acc_mes(n) + 2*sum(acc_mes(2:n-1)))*(h/2);
end

Diff1 = diff(acc_mes(1:end))/t(2);
Erreur = ((h^2)/12) * (Diff1(end) - Diff1(1))

%Mettre dans structure
Vitesse = struct("Y", Vitesse_mes, "X", t, "Erreur", Erreur);

%% Méthode Simpson
Hauteur_mes(1) = H_ini_NASA;
Position_mes(1) = t(1);
h = t(2);
for n = 3:2:length(Vitesse_mes)
    Hauteur_mes(((n-1)/2)+1,1) = Hauteur_mes(1) - (Vitesse_mes(1) + Vitesse_mes(n) + 4*sum(Vitesse_mes(2:2:(n-1))) + 2*sum(Vitesse_mes(3:2:(n-1)))) * (h/3);
    Position_mes(((n-1)/2)+1,1) = t(n);
end

Diff2 = diff(Diff1)/t(2);
Erreur = ((h^4)/180) * (Diff2(end) - Diff2(1))

%Mettre dans structure
Position = struct("Y", Hauteur_mes, "X", Position_mes, "Erreur", Erreur);

clear Vitesse_mes h acc_mes t Hauteur_mes Position_mes Erreur

%Affichage des approximations
figure
title("Vitesse")
scatter(Vitesse.X,Vitesse.Y)
grid on

figure
title("Position")
scatter(Position.X,Position.Y)
grid on

% À identifier
% RO_0
% h_s