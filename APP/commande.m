    function f = commande(t,z)
    
    %z = [v_ini, Gamma_ini, h_ini, s_ini, Theta_ini, q_ini]
    %

    load("variables.mat")

    func_v = z(1);
    func_Gamma = z(2);
    func_h = z(3);
    func_s = z(4);
    func_Theta = z(5);
    func_q = z(6);

    %Rayon
    func_r = R_mars + func_h; %Rayon
    func_g = (U_mars/(func_v*(func_r.^2)))';

    %Pdyn
    func_p = p0 * exp(-func_h/hs);
    func_P_dyn = (0.5) * func_p * func_v^2;

    %Gamma_ref
    func_Delta_V_Aero = V_fin(indice_gamma) - func_v;
    func_Gamma_ref = asind((0.5)*B*hs*((p_fin - func_p)./(log(1 + (func_Delta_V_Aero./func_v)))));

    %Alpha
    func_Alpha = func_Theta - func_Gamma;

    %Calcul de L_aero
    func_L_aero = func_P_dyn*S*C_Lalpha*func_Alpha;

    %Calcul de D_aero
    func_D_aero = func_P_dyn*S*C_D0;
    
    %Pour les theta commande
    temp1_1 = ((func_P_dyn*S*C_Lalpha*func_Gamma) / (func_v*m));
    temp2_1 = (((func_v./func_r)-func_g) * cosd(func_Gamma));
    temp3_1 = (K_p_trans*(func_Gamma_ref-func_Gamma));
    temp4_1 = (func_P_dyn*S*C_Lalpha) / (func_v*m);
    
    %Calcul de theta commande
    func_Theta_cmd = (temp1_1-temp2_1+temp3_1)./(temp4_1);

    %Ajustement de theta commande
    if func_Theta_cmd < -60
        func_Theta_cmd = -60;
    elseif func_Theta_cmd > 60
        func_Theta_cmd = 60;
    end

    %Gamma point
    if Asservissement == 1
        func_Gamma_point = (-temp1_1+temp2_1) + temp4_1*func_Theta_cmd;
    else
        func_Gamma_point = (1/func_v) * ((func_L_aero/m)*(((func_v^2)/func_r)-func_g)*cosd(func_Gamma));
    end

    %Pour les delta commande
    temp1_2 = (-1)*(((func_P_dyn*S*d*C_Malpha*func_Alpha)/J)+(func_P_dyn*S*d*C_Mq*func_q*(d/(2*J*func_v))));    
    temp2_2 = (K_p_rot) * (func_Theta_cmd - func_Theta);
            %pas sur de ça theta commande point c'est quoi 
    temp3_2 = (K_d_rot) * (func_Theta_cmd - func_q);
    temp4_2 = ((func_P_dyn*S*d*C_Mdelta)/J);
    
    %Calcul de delta commande
    func_Delta_cmd = (temp1_2 + temp2_2 + temp3_2)/temp4_2;

    %Gamma point
    if Asservissement == 1
        func_q_point = (-temp1_2) + temp4_2*func_Delta_cmd;
    else
        func_q_point = (1/J)*(func_P_dyn*S*d * ((C_Malpha*func_Alpha) + ((d*C_Mq*func_q)/(2*func_v))));
    end

    %Ressortir toute les valeurs qui sont nécessaire
    f(1) = ((-1)*(func_D_aero/m)) - (func_g*sind(func_Gamma));
    f(2) = func_Gamma_point;
    f(3) = func_v * sind(func_Gamma);
    f(4) = (func_v/func_r) * cosd(func_Gamma);
    f(5) = func_q;
    f(6) = func_q_point;

    f = f(:);

    %Gamma
    %Gamma_ref
    %v
    %

end

