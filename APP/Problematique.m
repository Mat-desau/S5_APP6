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
C_Mdelta = 0.10;
B = (C_D0*S)/m;

%% Conditions initiales
V_ini = 6100; %m/s
Gamma_ini = -20.5; %degré
h_ini = 120000; %m
s_ini = 0.0; %degré
Theta_ini = -80; %deg
q_ini = 0.0; %degré/s
indice_gamma = 1;   %1 = 250m/s     2 = 300m/s
Asservissement = 1; %1 = on         0 = off
z0 = [V_ini, deg2rad(Gamma_ini), h_ini, deg2rad(s_ini), deg2rad(Theta_ini), deg2rad(q_ini), 0];

%% Conditions finales désirer
V_fin1 = 250; %m/s
V_fin2 = 300; %m/s
V_fin = [250 300]; %m/s
h_fin = 10000; %m

%% Conditions initiales NASA
Gamma_ini_NASA = -90; %degrée
s_ini_NASA = 0.0; %degrée
Theta_ini_NASA = -90; %degrée
q_ini_NASA = 0.0; %degrée

%% Limites Strucutrel
P_dyn_max = 9500;  %N/m^2
Delta_t_lim = 40; %sec
D_aero_max = 2650; %N
Delta_Cmd_max = 60; %degrée

%% Load les valeurs
load("Accelero_Data_from_NASA.mat")
t_28 = t(1:2:end);

%% Accélération
%Mettre dans structure
Accelertion = struct("X", t, "Y", acc_mes);

%% Méthode de trapèze
%On commence a 6100
vit_mes(1,1) = V_ini;
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
h_mes(1) = h_ini;
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

%Transposition des valeurs
m_cal = A(2);
b_cal = A(1);

%Voir calculs à la main sur comment on a trouver

%Identification des paramètres manquant
p0 = (m.*exp(b_cal))./(S.*C_D0);
hs = -1./m_cal;

%Approximation de l'accélération
acc_appro = (0.5 * p0 * S * C_D0 .* exp(-h_mes/hs) .* vit_mes(1:2:end).^2)/m;

%Erreur
Y_appro = log((2*acc_appro)./(vit_mes(1:2:end)).^2);
N = length(acc_mes);
y_ = (1/N) * sum(Y);
R_2_iden = (sum((Y_appro - y_).^2) ./ sum((Y - y_) .^2)); 
clear y_ N

%Calcul des erreurs
N = length(Y);
y_ = (1/N) * sum(Y);
R_2_acc = (sum((Y_appro-y_).^2) ./ sum((Y - y_) .^2))
clear y_ N

N = length(acc_appro);
RMS_abs_acc = sqrt((1/N) * sum((acc_appro - acc_mes(1:2:end)).^2));
RMS_rel_acc = sqrt((1/N) * sum(((acc_appro - acc_mes(1:2:end))./acc_mes(1:2:end)).^2));
clear N

%% Conception des asservissements RAA
%Valeurs
p_ini = p0 * exp(-h_ini/hs);
p = p0 * exp(-h_mes/hs);

%Calcul de la formule
vit_appro = V_ini * exp(0.5*B*hs*((p-p_ini)/sind(Gamma_ini_NASA)));

%Nouveau calcul de l'accélération par rapport a vitesse approximé
acc_appro_vit_appro = (0.5 * p0 * S * C_D0 .* exp(-h_mes/hs) .* vit_appro.^2)/m;

%% Limites structurel
p = p0 * exp(-h_ini/hs);
p_fin = p0 * exp(-h_fin/hs);

%Delta_V_Aero(1) = V_fin1 - V_ini
%Delta_V_Aero(2) = V_fin2 - V_ini 
Delta_V_Aero(1) = V_fin(1) - sqrt((V_ini^2)+(2*U_mars*((1/(R_mars+h_fin))-(1/(R_mars+h_ini)))));
Delta_V_Aero(2) = V_fin(2) - sqrt((V_ini^2)+(2*U_mars*((1/(R_mars+h_fin))-(1/(R_mars+h_ini)))));

%Calcul de Gamma
Gamma_ref(1) = asind((0.5)*B*hs*((p_fin - p)./(log(1 + (Delta_V_Aero(1)./V_ini))))); %250m/s
Gamma_ref(2) = asind((0.5)*B*hs*((p_fin - p)./(log(1 + (Delta_V_Aero(2)./V_ini))))); %300m/s

%RAA pour trouver vitesse
h_array = [h_ini:-10:h_fin]';
Distance = h_ini - h_array;
p_ini = p0 * exp(-h_ini/hs);
p_array = p0 * exp(-h_array/hs);

vit_cal_RAA(:,1) = V_ini * exp(0.5*B*hs*((p_array-p_ini)/sind(Gamma_ref(1)))); %250m/s
vit_cal_RAA(:,2) = V_ini * exp(0.5*B*hs*((p_array-p_ini)/sind(Gamma_ref(2)))); %250m/s

P_dyn_cal(:,1) = (0.5)*p_array.*(vit_cal_RAA(:,1).^2);
P_dyn_cal(:,2) = (0.5)*p_array.*(vit_cal_RAA(:,2).^2);

D_aero_cal(:,1) = P_dyn_cal(:,1)*S*C_D0;
D_aero_cal(:,2) = P_dyn_cal(:,2)*S*C_D0;

%% Calcul avec Newton-Raphson Limites Structurel
%Pour trouver ça on a mit 
h_dep(1) = 34000;
h_dep(2) = 20000;

for n = 1:2
    x = h_dep(n);

    %Calcul de la dériver
    P_vrai = p0*exp(-x/hs);
    P_prime = (-p0/hs)*exp(-x/hs);
    
    V_vrai = V_ini*exp(0.5*B*hs*((P_vrai-p_ini)/sind(Gamma_ref(indice_gamma))));
    V_prime = (((0.5)*B*hs*V_ini)/sind(Gamma_ref(indice_gamma))) * P_prime * exp(0.5*B*hs*((P_vrai-p_ini)/sind(Gamma_ref(indice_gamma)))); %avec 250 km/h
    
    % On enleve D_aero_max pour crée des nouveau zeros
    D_aero_vrai = ((0.5) * P_vrai * (V_vrai.^2) * S * C_D0) - D_aero_max;
    D_aero_prime = ((0.5) * S * C_D0 * ((P_prime * (V_vrai.^2))+(2*P_vrai*V_vrai*V_prime)));
    
    %Pour calcul Newton-Rapson
    tol = 1*10^(-10);
    iter = 0;
    Max_iter = 1000;
    
    %Boucle qui cherche
    while (iter < Max_iter) && (abs(D_aero_vrai) > tol)
        x = x - D_aero_vrai/D_aero_prime;
        
        P_vrai = p0*exp(-x/hs);
        P_prime = (-p0/hs)*exp(-x/hs);
    
        V_vrai = V_ini*exp(0.5*B*hs*((P_vrai-p_ini)/sind(Gamma_ref(indice_gamma))));
        V_prime = (((0.5)*B*hs*V_ini)/sind(Gamma_ref(indice_gamma))) * P_prime * exp(0.5*B*hs*((P_vrai-p_ini)/sind(Gamma_ref(indice_gamma)))); %avec 250 km/h
        
        % On enleve D_aero_max pour faire que les "zeros" sont a D_aero_max
        D_aero_vrai = ((0.5) * P_vrai * (V_vrai.^2) * S * C_D0) - D_aero_max;
        D_aero_prime = ((0.5) * S * C_D0 * ((P_prime * (V_vrai.^2))+(2*P_vrai*V_vrai*V_prime)));
    
        iter = iter + 1;
    end
    
    %Iter - Tol - D_aero - Vitesse - Position
    Val(n,1) = iter;
    Val(n,2) = abs(D_aero_vrai);
    Val(n,3) = D_aero_vrai;
    Val(n,4) = V_vrai;
    Val(n,5) = x;
end

%Avec ces points on peut trouver les autres données
h_min = Val(2,5);
h_max = Val(1,5);

v_min = Val(2,4);
v_max = Val(1,4);

vit_moy = (v_max + v_min) / 2;

Delta_t = abs((abs(h_max-h_min)/(vit_moy*sind(Gamma_ref(indice_gamma)))));


%% Commande dynamique
%Pour Kp
Tau_trans = 0.25;
K_p_trans = 1/Tau_trans;

Zeta_rot = 0.7;
Omega_rot = 20; %rad/s

K_p_rot = Omega_rot^2;
K_d_rot = 2*Zeta_rot*Omega_rot;

%Sortir les variables nécessaires pour le code
save variables.mat p0 R_mars U_mars hs V_fin p_fin indice_gamma B S C_Lalpha m K_p_trans C_Malpha d J C_Mq K_p_rot K_d_rot C_Mdelta C_D0 Asservissement h_fin

tspan = [0 300];
reltol2 = 1e-10;
options = odeset('abstol', 1e-06, 'reltol', reltol2);
[func_t, func_z] = ode45('commande', tspan, z0, options);

index = find(func_z(:,3) <= h_fin);


%Valeur nécessaire après la sortie de Valid
%P_dyn
P_dyn_Valid = (0.5) * (p0*exp(-func_z(:,3)./hs)) .* (func_z(:,1).^2);
P_dyn_max_Valid = max(P_dyn_Valid);
%Gamma ref
func_p_Valid = p0 * exp((-func_z(:,3))./hs);
Delta_V_Aero_Valid = V_fin(indice_gamma) - sqrt((func_z(:,1).^2)+(2*U_mars.*((1/(R_mars+h_fin))-(1./(R_mars+func_z(:,3))))));
Gamma_ref_Valid = asin((0.5)*B*hs*((p_fin - func_p_Valid)./(log(1 + (Delta_V_Aero_Valid./func_z(:,1))))));
%D_aero
D_aero_Valid = P_dyn_Valid*(S*C_D0);
%Delta_t
Delta_t_Valid = func_z(:,7);
%Alpha
Alpha_Valid = func_z(:,5) - func_z(:,2);


%% Texte
% disp("---------------- INTEGRATION ----------------");
% texte = ["Erreur Intégration Vitesse : ", Vitesse.Erreur];
% disp(texte);
% texte = ["Erreur Intégration Position : ", Position.Erreur];
% disp(texte);
% 
% disp("---------------- ACCELERATION ----------------");
% texte = ["Valeur de RMS_abs accélération : ", RMS_abs_acc];
% disp(texte);
% texte = ["Valeur de RMS_rel accélération : ", RMS_rel_acc];
% disp(texte);
% texte = ["Valeur de R^2 accélération : ", R_2_acc];
% disp(texte);
% 
% disp("---------------- IDENTIFICATION ----------------");
% texte = ["m : ", m_cal];
% disp(texte);
% texte = ["b : ", b_cal];
% disp(texte);
% texte = ["P_0 : ", p0];
% disp(texte);
% texte = ["h_s : ", hs];
% disp(texte);
% texte = ["Valeur de R^2 Identification : ", R_2_iden];
% disp(texte);
% 
% disp("---------------- LIMITES STRUCTUREL ----------------");
% disp("Gamma")
% texte = ["Gamma avec 250m/s : ", Gamma_ref(1)];
% disp(texte);
% texte = ["Gamma avec 300m/s : ", Gamma_ref(2)];
% disp(texte);
% disp("P_dyn");
% texte = ["P_dyn max avec 250m/s : ", max(P_dyn_cal(:,1))];
% disp(texte);
% texte = ["P_dyn max avec 300m/s : ", max(P_dyn_cal(:,2))];
% disp(texte);
% 
% disp("Newton-Raphson pour indice gamma = (1 = 250m/s)(2 = 300m/s)");
% disp(indice_gamma);
% 
% texte = ["h min avec indice_gamma : ", h_min];
% disp(texte);
% texte = ["v min avec indice_gamma : ", v_min];
% disp(texte);
% texte = ["h_depart avec indice_gamma : ", h_dep(2)];
% disp(texte);
% texte = ["Iteration avec indice_gamma : ", Val(2,1)];
% disp(texte);
% 
% texte = ["h max avec indice_gamma : ", h_max];
% disp(texte);
% texte = ["v max avec indice_gamma : ", v_max];
% disp(texte);
% texte = ["h_depart avec indice_gamma : ", h_dep(1)];
% disp(texte);
% texte = ["Iteration avec indice_gamma : ", Val(1,1)];
% disp(texte);
% 
% texte = ["Delta_t avec indice_gamma : ", Delta_t];
% disp(texte);
% 
% disp("---------------- VALID ----------------");
% disp("Asservissement -> 1=on 2=off")
% Asservissement
% disp("V_ini -> 1=250 2=300 ")
% indice_gamma
% 
% texte = ["V_finale : ", func_z(index(1),1)];
% disp(texte);
% texte = ["Delta_t : ", max(func_z(:,7))];
% disp(texte);
% texte = ["P_dyn_max : ", P_dyn_max_Valid];
% disp(texte);


%% Graphiques
% -------------------- Lissage -------------------------
%Graphique Accélération approximée par lissage
% figure
% hold on
% scatter(t, acc_mes, "black")
% plot(t_28, acc_appro, "red")
% grid on
% title("Accélération approximée par lissage")
% legend(["Mesurer" "Approximer"])

%Graphique Vitesse obtenue par intégration
% figure
% hold on
% scatter(Vitesse.X, Vitesse.Y, "black")
% grid on
% title("Vitesse obtenue par intégration")
% legend(["Mesurer"])

%Graphique Altitude obtenue par intégration
% figure
% hold on
% scatter(Position.X, Position.Y, "black")
% grid on
% title("Altitude obtenue par intégration")
% legend(["Mesurer"])

% -------------------- RAA -------------------------

%Graphique Vitesse calculer par RAA
% figure
% hold on
% scatter(h_mes, vit_mes(1:2:end), "blue", 'filled')
% plot(h_mes, vit_appro, "red")
% grid on
% xlabel("Hauteur (m)")
% ylabel("Vitesse (m/s)")
% title("Vitesse calculé par RAA")
% legend(["Mesurer" "Approximer"])

%Graphique Accélération calculé par RAA
% figure
% hold on
% scatter(h_mes, acc_mes(1:2:end), "blue", 'filled')
% plot(h_mes, acc_appro_vit_appro, "red")
% grid on
% xlabel("Hauteur (m)")
% ylabel("Accélération (m/s^2)")
% title("Accélération calculé par RAA")
% legend(["Mesurer" "Approximer"])

% -------------------- Limitation Structurel -------------------------

%Graphique Vitesse avec les limitations structurel
% figure
% hold on
% plot(Distance, vit_cal_RAA)
% grid on
% title("Vitesse calculé par RAA")
% legend(["250m/s" "300m/s"])

%Graphique P_dyn avec les limitations structurel
% figure
% hold on
% plot(Distance, P_dyn_cal(:,1), "red")
% plot(Distance, P_dyn_cal(:,2), "blue")
% plot(Distance, P_dyn_max*ones(length(Distance),1), "black")
% grid on
% title("P_d_y_n calculé par RAA")
% legend(["250m/s" "300m/s" "P_d_y_n max"])

%Graphique D_aero avec les limitations structurel
% figure
% hold on
% plot(Distance, D_aero_cal(:,1), "red")
% plot(Distance, D_aero_cal(:,2), "blue")
% plot(Distance, D_aero_max*ones(length(Distance),1), "black")
% grid on
% title("D_a_e_r_o calculé par RAA")
% legend(["250m/s" "300m/s" "D_a_e_r_o max"])

% -------------------- Validation -------------------------

%Graphique des asservissements
% figure;
% hold on
% grid on
% plot((func_z(1:index(1),3)/1000), func_z(1:index(1),1));
% xlabel("Hauteur (km)");
% ylabel("Vitesse (m/s)");
% legend(["Vitesse"])
% title("Vitesse en fonction de la hauteur");
% 
% figure;
% hold on
% grid on
% plot((func_t(1:index(1))), func_z(1:index(1),1));
% xlabel("Temps (s)");
% ylabel("Vitesse (m/s)");
% legend(["Vitesse"])
% title("Vitesse en fonction du temps");
% 
% if Asservissement == 1
    % figure;
    % hold on
    % grid on
    % plot(func_t(1:index(1)), rad2deg(func_z(1:index(1),2)));
    % plot(func_t(1:index(1)), rad2deg(Gamma_ref_Valid(1:index(1))))
    % % plot(func_t(1:index(1)), Gamma_ref(indice_gamma)*ones(index(1),1), "black")
    % % legend(["\gamma" "\gamma_r_e_f" "\gamma_r_e_f initial"])
    % legend(["\gamma" "\gamma_r_e_f"])
    % xlabel("Temps (s)");
    % ylabel("Angle (deg)");
    % ylim([-25 -10])
    % title("Gamma en fonction du temps");
% 
%     figure;
%     hold on
%     grid on
%     plot(func_z(1:index(1),3)/1000, rad2deg(func_z(1:index(1),2)));
%     legend(["\gamma"])
%     xlabel("Hauteur (km)");
%     ylabel("Angle (deg)");
%     title("Gamma en fonction de hauteur");
% else
%     figure;
%     hold on
%     grid on
%     plot(func_t(1:index(1)), rad2deg(func_z(1:index(1),2)));
%     legend(["\gamma"]);
%     xlabel("Temps (s)");
%     ylabel("Angle (deg)");
%     title("Gamma en fonction du temps");
% end
% 
% figure;
% hold on
% grid on
% plot(func_t(1:index(1)), func_z(1:index(1),3)/1000);
% xlabel("Temps (s)");
% ylabel("Hauteur (km)");
% legend(["Hauteur"]);
% title("Hauteur en fonction du temps");
% 
% figure;
% hold on
% grid on
% plot(func_t(1:index(1)), rad2deg(func_z(1:index(1),5)));
% plot(func_t(1:index(1)), rad2deg(Alpha_Valid(1:index(1))));
% xlabel("Temps (s)");
% ylabel("Angle (deg)");
% legend(["\Theta" "\alpha"]);
% title("Theta et Alpha");
% 
if Asservissement == 1
    figure;
    hold on
    grid on
    plot(func_t(1:index(1)), rad2deg(func_z(1:index(1),6)));
    xlabel("Temps (s)");
    ylim([-5 5])
    ylabel("Angle (deg)");
    legend(["q"]);
    title("q");
else
    figure;
    hold on
    grid on
    plot(func_t(1:index(1)), rad2deg(func_z(1:index(1),6)));
    xlabel("Temps (s)");
    ylabel("Angle (deg)");
    ylim([-5 5])
    legend(["q"]);
    title("q");
end
% 
% figure;
% hold on
% grid on
% plot(func_t(1:index(1)), P_dyn_Valid(1:index(1)));
% plot(func_t(1:index(1)), D_aero_Valid(1:index(1)));
% xlabel("Temps (s)");
% ylabel("Amplitude");
% legend(["P_d_y_n" "D_a_e_r_o"]);
% title("D_a_e_r_o et P_d_y_n");
% 
% figure;
% hold on;
% grid on;
% plot(func_t(1:index(1)), Delta_t_Valid(1:index(1)));
% ylabel("Temps dépassement (s)");
% xlabel("Temps chute (s)")
% title("Delta temps")
% legend(["\Delta t"])


disp("Hello World")