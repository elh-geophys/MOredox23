% EL
% Original: Aug 2022
% Updated: July 2024
%      - updated FeO tracking and dIW
%      - updated with new compositions
%
% Function for getting metal rain ratio with varying depth and
% efficiencies.  Use for mass sampling.
%
% INPUTS:
%
%           j indexing is function of time
%           i indexing is function of mantle depth
%
%   P_cmb         [Pa] pressure at CMB during accretion, i length array
%   P_min         [Pa] pressure min for MO depth with only impactor melt
%   dP            [Pa] increments of P
%   r_0           [] initial Fe3/Fe ratio at t=0 for Earth
%   Accr_model    [] earth mass growth model; j length array
%   dMp_temp      [kg] the potential mass equilibrated during time interval dt; 
%                       j-1 length array
%   compSheet_earth   [] sheet to use for Earth composition, j length array
%   FeO_imp       [] the FeO wt % of impactor, j length array
%   r_imp         [] Fe3+/sumFe ratio for impactor, j length array
%   rho_m         [kg/m^3] density of mantle, j length array
%   M_m           [kg] mass of mantle, j length array
%   M_m_imp       [kg] mass of impactor silicate, j-1 length array
%   R_E_post      [m] radius of Earth post impact, j-1 array
%   R_c_post      [m] radius of Earth's core post impact, j-1 array
%   Tp_const      [K] constant temperature Tp value


function [r_m_Dt, eff_out, P_mo_out] = getRainRatio_constT(P_cmb, P_min, dP, r_0, Accr_model, ...
        dMp_temp, compSheet_earth, FeO_imp, r_imp, rho_m, rho_c_post, M_m, M_m_imp, R_E_post, R_c_post, Tp_const)
    
    P_cmb_max = round(P_cmb(end),-9);
    P_length_max = P_cmb_max/dP+1;             %to go by dP and +1 for end points :)

    r_m_Dt = zeros(1,length(Accr_model));      % final Fe3+/sumFe of mantle over evolution time
    r_m_Dt(1) = r_0;
    r_m = zeros(P_length_max+1,1);        %temp Fe3+/sumFe for each iteration of fall time
    
    eff_out = zeros(1, length(Accr_model));      % to record random efficiencies
    P_mo_out = zeros(1, length(Accr_model));     % to record random depth of melting
    
    PV_data = readmatrix('/db/PVcalc.xlsx');
    Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');
    CompEarth_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_earth);
    MolW_byM_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Emantle', 'Range', 'B2:B13');
    MolW_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Emantle', 'Range', 'C2:C13');
    
    FeO_E = CompEarth_data(7,:);

    for j = 1:length(Accr_model)-1      % earth evolution time

        if round(dMp_temp(j),3) == 0    % NOT AN IMPACT
            r_m_Dt(j+1) = r_m_Dt(j);    % ratio stays the same at this time step

        else                            % IMPACT, VARIANT MO DEPTH AND EFFICIENCY
            
            % CHOOSE MAGMA OCEAN DEPTH
            eff = rand(1);              %[] randomize efficiency factor for droplet equilibrium
            eff_out(j) = eff;
            d_mo_factor = rand(1);      %[] randomize depth of melting
            P_mo = d_mo_factor*(P_cmb(j)-P_min(j))+P_min(j);        % P_mo between P_min and P_cmb
            P_mo_out(j) = P_mo;
            
            P_max = min(round(P_mo,-9), P_cmb(j)-mod(P_cmb(j),0.5e9));         %set up P, with rare case where P_mo rounds up above P_cmb
            P = (0:dP:P_max)';              
            if size(P)==1
                P = [0, dP];                %for some situations where P very small
            end
           
            % DETERMINE MAGMA OCEAN GEOTHERM
            % choice where Tp is constant          
            Tp = Tp_const;
            Tad = getMOAdiabat(Tp,P, Adiabat_data);
            
            % DETERMINE PV INTEGRATION AS INT(PV)/RT
            PV = calcPV(Tad,P,PV_data);

            % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
            [r_eq,~] = calcFeRatio(Tad, P, PV, CompEarth_data(3:end,j+1), MolW_data, MolW_byM_data);
            
            % METAL RAIN CALCULATION FOR THIS TIME STEP
            R = linspace(R_c_post(j),R_E_post(j),1000);
            P_check = get2LayerP(rho_m(j+1), rho_c_post(j), R_E_post(j), R_c_post(j), R);
            if P_mo > P_check(1)
                P_mo = P_check(1);   %the rare floating point issue, when d_mo_factor = 1 and P_mo isn't exactly P_check(1) as it should be
            end
            R_mo = interp1(P_check, R, P_mo);
            M_mo = rho_m(j+1) * 4/3*pi*(R_E_post(j)^3 - R_mo^3);        %approximate mass of spherical shell MO
        
            % initial Fe ratio of mantle after mixing previous silicate with impactor silicate
            r_m(1) = ((M_mo-M_m_imp(j))*FeO_E(j)*r_m_Dt(j)+M_m_imp(j)*FeO_imp(j)*r_imp(j))/((M_mo-M_m_imp(j))*FeO_E(j) + M_m_imp(j)*FeO_imp(j));
            FeO_mo = ((M_mo-M_m_imp(j))*FeO_E(j) + M_m_imp(j)*FeO_imp(j))/M_mo;
                
            dMp = min(dMp_temp(j)*eff, M_mo);            %silicate mass eq'd with given efficiency, if eq'd mass > spherical shell mass, take all MO mass eq'd

            for i = 1:length(P)          % droplets falling through mantle
                if P(i) < P_mo
                    % ratio changes as the droplets fall through mantle layers until base of MO
                    r_m(i+1) = (dMp*r_eq(i) + (M_mo-dMp)*r_m(i))/(M_mo);
                    idx = i+1;
                end
            end

            % final r from mixing redox'd MO with whole mantle
            r_m_Dt(j+1) = (M_mo*FeO_mo*r_m(idx)+(M_m(j+1)-M_mo)*FeO_E(j)*r_m_Dt(j))/(M_mo*FeO_mo + (M_m(j+1)-M_mo)*FeO_E(j));
                
            r_m = r_m*0;                % reset r_m
            
        end
    end  
end

