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
%Mettre dans structure
Accelertion = struct("X", t, "Y", acc_mes);

%% Méthode de trapèze
%On commence a 6100
vit_mes(1,1) = V_ini_NASA;
Delta_x = t(2);

%Faire la fonction
for n = 2:length(acc_mes)
    vit_mes(n,1) = vit_mes(1,1) - (acc_mes(1) + acc_mes(n) + 2*sum(acc_mes(2:n-1)))*(Delta_x/2);
end

%Differentiel pour erreur
Diff1(1) = (acc_mes(2)-acc_mes(1))/Delta_x;
Diff1(2) = (acc_mes(end)-acc_mes(end-1))/Delta_x;

Erreur = ((Delta_x^2)/12) * (Diff1(2) - Diff1(1));

%Mettre dans structure
Vitesse = struct("Y", vit_mes, "X", t, "Erreur", Erreur);

%% Méthode Simpson
%On commence a une hauteur initiale
h_mes(1) = H_ini_NASA;
t_h_mes(1) = t(1);
Delta_x = t(2);

%Faire la fonction
for n = 3:2:length(vit_mes)
    h_mes(((n-1)/2)+1,1) = h_mes(1) - (vit_mes(1) + vit_mes(n) + 4*sum(vit_mes(2:2:(n-1))) + 2*sum(vit_mes(3:2:(n-1)))) * (Delta_x/3);
    t_h_mes(((n-1)/2)+1,1) = t(n);
end

%Différentiel pour erreur
Diff3(1) =  (vit_mes(4) - 3*vit_mes(3) + 3*vit_mes(2) - vit_mes(1))/ (Delta_x^3);
Diff3(2) =  (vit_mes(end) - 3*vit_mes(end-1) + 3*vit_mes(end-2) - vit_mes(end-3))/ (Delta_x^3);

Erreur = ((Delta_x^4)/180) * (Diff3(2) - Diff3(1));

%Mettre dans structure
Position = struct("Y", h_mes, "X", t_h_mes, "Erreur", Erreur);

%% Identification p0 et hs
%Moindre Carrer
C = [ones(length(h_mes),1), h_mes];
Y = log((2*acc_mes(1:2:end))./(vit_mes(1:2:end)).^2);
A = pinv(C)*Y;

%Voir calculs à la main sur comment on a trouver

%Identification des paramètres manquant
p0 = (m.*exp(A(1)))./(S.*C_D0);
hs = -1./A(2);

%Approximation de l'accélération
acc_appro = (0.5 * p0 * S * C_D0 .* exp(-h_mes/hs) .* vit_mes(1:2:end).^2)/m;

%Graphique Acceleration approx vs mesurer
% figure
% hold on
% scatter(t, acc_mes, "red")
% plot(t_h_mes, acc_appro, "blue")
% grid on
% title("Accélération")
% legend(["Mesurer" "Approximer"])

%Calcul des erreurs
N = length(acc_mes);
y_ = (1/N) * sum(acc_mes);
R_2 = (sum((acc_appro-y_).^2) ./ sum((acc_mes - y_) .^2));
clear y_ N

N = length(acc_appro);
RMS_abs = sqrt((1/N) * sum((acc_appro - acc_mes(1:2:end)).^2));
RMS_rel = sqrt((1/N) * sum(((acc_appro - acc_mes(1:2:end))./acc_mes(1:2:end)).^2));
clear N

%% Conception des asservissements
%Valeurs
B = (C_D0*S)/m;
p_ini = p0 * exp(-h_ini/hs)
p = p0 * exp(-h_mes/hs)

%Calcul de la formule
vit_appro = V_ini * exp(0.5*B*hs*((p-p_ini)/sin(Theta_ini_NASA)))

%Graphique Vitesse approx vs mesurer
% figure
% hold on
% scatter(t, vit_mes, "red")
% plot(t(1:2:end), vit_appro, "blue")
% grid on
% title("Vitesse")
% legend(["Mesurer" "Approximer"])



