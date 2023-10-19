% EL
% Original: July 2022
% Updated: 2023-09-26
%
% Iron droplet equilibrium through MO during accretion with given efficiency
%
% Change the model type, efficiency, depth of MO, comp sheets, Tp type,
% file writes
%
% This script can be used to test constant efficiencies (eff) and MO depths
% (d_mo_factor)

clear;

% PARAMETERS TO CHANGE
model = 4;                      %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
r_0 = 0.004;                    %initial Fe3+/sumFe for Earth
eff = 0.5;                        %[] efficiency factor for droplet equilibrium, 1 = 100%
d_mo_factor = 0.5;                %[] fraction of Pcmb pressure for MO base pressure
compSheet_early = 'EarthEarly';           %sheet in MoleWeights.xlsx to use for composition
compSheet_late = 'EarthLate';
Tp_type = 'constant';           %chooose method to calculate Tp, either 'Pmo', 'U2Q', or 'constant'
    T0 = 1613;                      %[K] initial Tp if using U2Q method
    epsilon = 0.2;                  %energy contribution factor for U
    Tp_const = 3500;                %[K] Tp for 'constant' method
dP = 0.5e9;                     %[Pa] increments of P for layer to do metal rain calculation

sheetOut = 'H04';               %sheet name to record data
fileOut = 'Rain.xlsx';          %file name to record data
write = 0;                      %1 or 0, to write to file

% ---------------------------------------------------------------------- %

% READ DATA SHEETS
PV_data = readmatrix('/db/PVcalc.xlsx');
Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');
CompEarly_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet_early);
CompLate_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet_late);

% CONSTANTS
Cp = 1e3;               %[J/kgK] specific heat 
M_E_0 = 5.97e24;        %[kg] present day mass of Earth
M_c_0 = 1.88e24;        %[kg] present day mass of core
rho_imp = 5000;         %[kg/m^3] approximation based on weighted average (0.68Si + 0.32Fe)

% SET UP ACCRETION MODEL
[t,Accr_model] = getAccrModel(model);
M_E = M_E_0 * Accr_model;
M_imp = M_E_0 * diff(Accr_model);

% note: during impact "n", earth mass is "n" pre-impact and "n+1" post-impact

% assume core and mantle take up proportional mass of Earth
M_c = M_E * M_c_0/M_E_0;    %core
M_m = M_E - M_c;            %mantle

% Rubie+2011, approximations for mantle and core densities
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

V_c = M_c./rho_c;                           %volume of core
V_m = M_m./rho_m;                           %volume of mantle
R_E = (3/(4*pi)*(V_m + V_c)).^(1/3);        %radius of Earth
R_imp = (3/(4*pi)*M_imp/rho_imp).^(1/3);    %radius of impactor 

Fe_ratio = 0.321;               % [] weight % of iron on Earth/impactor
Si_ratio = 1-Fe_ratio;          % [] weight % of silicate on Earth/impactor
M_c_imp = M_imp*Fe_ratio;       % [kg] approximate proportion metal mass of impactor
M_m_imp = M_imp*Si_ratio;       % [kg] approximate proportion silicate mass of impactor

%estimated Fe3+/sumFe for impactor, based on Rubie+2011, Supp. Table 3a
%determined endpoints by using test_getSingleFeRatio_H22.m
r_imp = zeros(1, length(M_imp));
for i = 1:length(M_imp)
    if Accr_model(i+1) <= 0.6
        % so GI of N21 model will have the higher value
        r_imp(i) = 0.004;                   
    else
        r_imp(i) = 0.0122;
    end
end

% Use half impactor core mass between pre- and post-impact to determine R_c at impact
% Use entire mantle mass for chemical mixing post-impact
M_c_post = M_c(1:end-1)+M_c_imp/2;                      %core mass with half of the impactor 
rho_c_post = (rho_c(1:end-1)+rho_c(2:end))/2;           %average core density between pre and post impact
R_c_post = (3/(4*pi)*(M_c_post./rho_c_post)).^(1/3);    %radius of core post impact, Earth + half of impactor

% radius of Earth w/ all impactor Si (M_m post-impact) and 1/2 impactor Fe
R_E_post = (3/(4*pi)*(M_m(2:end)./rho_m(2:end) + M_c_post./rho_c_post)).^(1/3);

% radius of Earth w/ no impactor Si (M_m pre-impact), but 1/2 impactor Fe;
% upper bound for MO radius with all impactor as melt
R_E_post_noimp = (3/(4*pi)*(M_m(1:end-1)./rho_m(2:end) + M_c_post./rho_c_post)).^(1/3);

%T&S Geodyanm Eqn 2.73, Pressure as a function radius from center
% determine pressure at CMB
P_cmb = zeros(1,length(M_imp));
for i = 1:length(M_imp)
    P_cmb(i) = get2LayerP(rho_m(i+1), rho_c_post(i), R_E_post(i), R_c_post(i), R_c_post(i));
end
P_cmb(isnan(P_cmb))=0;          %for some accretion models, where P=0 at t=0
P_cmb_max = round(P_cmb(end),-9);
P_length_max = P_cmb_max/dP+1;             %to go by dP and +1 for end points :)

% determine the minimum pressure without impactor melt
P_min = zeros(1,length(M_imp));
for i = 1:length(M_imp)
    P_min(i) = get2LayerP(rho_m(i+1), rho_c_post(i), R_E_post(i), R_c_post(i), R_E_post_noimp(i));
end
P_min(isnan(P_min))=0;          %for some accretion models, where P=0 at t=0

[M_eq] = calcEqSi(M_c_imp, 'sph');
dMp_temp = M_eq*eff;

r_m_Dt = zeros(1,length(t));        %final Fe3+/sumFe of mantle over evolution time
r_m_Dt(1) = r_0;
r_m = zeros(P_length_max+1,1);        %temp Fe3+/sumFe for each iteration of fall time (+1 b/c of initial mixing of mantle & impactor)
r_m(1) = r_0;

for j = 1:length(t)-1                 % j index for evolution time, where j=impact point, j+1=after impact
    if M_imp(j) > 0
        
        % CHOOSE MAGMA OCEAN DEPTH
        P_mo = d_mo_factor*(P_cmb(j)-P_min(j))+P_min(j);
        
        P_max = round(P_mo,-9);         %set up P
        P = (0:dP:P_max)';              
        if size(P)==1
            P = [0, dP];                %for some accretion models where early t, P very small
        end
        
        % DETERMINE MAGMA OCEAN GEOTHERM
        Tp_mo_base = getMOTp(P_mo/1e9, Adiabat_data);
        
        switch Tp_type
            case 'Pmo'
                % choice where Tp is based on full MO being at above liquidus
                Tp = Tp_mo_base;
                
            case 'U2Q'
                % choice where Tp is based on increasing temperature from T0 due to impact GPE
                % if Tp turns out to be > Tp_mo_base, then take it all as molten 
                Tf_test = calcTforU2Q(T0,Cp,M_E(j),M_c(j),M_imp(j),R_E(j),R_imp(j),Si_ratio,epsilon);

                if Tf_test > Tp_mo_base
                    Tp = Tp_mo_base;
                else
                    Tp = Tf_test;
                end
                
            case 'constant'
                % choice where Tp is a constant value through all accretion
                Tp = Tp_const;
                
            otherwise
                disp('Could not calculate Tp based on input Tp_type')
                Tp = NaN;
                break               % ends the if statement
        end        
        
        Tad = getMOAdiabat(Tp,P, Adiabat_data);
        
        % DETERMINE PV INTEGRATION AS INT(PV)/RT
        PV = calcPV(Tad,P,PV_data);
        
        % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
        % Use post-impact comp and dIW
        [r_eq,dIW] = calcFeRatio(Tad, P, Accr_model(j+1), PV, CompEarly_data, CompLate_data);
        disp([num2str(round(Tp)), 'K and ', num2str(dIW)])
        
        % METAL RAIN CALCULATION FOR THIS TIME STEP
        R = linspace(R_c_post(j),R_E_post(j),1000);
        P_check = get2LayerP(rho_m(j+1), rho_c_post(j), R_E_post(j), R_c_post(j), R);
        R_mo = interp1(P_check, R, P_mo);
        M_mo = rho_m(j+1) * 4/3*pi*(R_E_post(j)^3 - R_mo^3);        %approximate mass of spherical shell MO
            
        % initial Fe ratio of mantle after mixing previous silicate with impactor silicate
        r_m(1) = ((M_mo-M_m_imp(j))*r_m_Dt(j)+M_m_imp(j)*r_imp(j))/(M_mo);
        
        dMp = min(dMp_temp(j), M_mo);       %take minimum between equilibrated mass & magma ocean mass
        
        for i = 1:length(P)                 % droplets falling through mantle time, r_m(1) = previous iteration's r_m(end)
            if P(i) < P_mo
                if isfinite(r_eq(i))        % only do mixing above T40, in the case of some U2Q method impacts
                    % ratio changes as the droplets fall through mantle layers until base of MO
                    % note r_m(i+1) is at P(i) location
                    r_m(i+1) = (dMp*r_eq(i) + (M_mo-dMp)*r_m(i))/(M_mo);
                    idx = i+1;
                end
            end
        end
        
        figure(4);
        hold on
        plot(r_eq, 'k--')
        plot(r_m, 'm')
        
        r_m_Dt(j+1) = (M_mo*r_m(idx) + (M_m(j+1)-M_mo)*r_m_Dt(j))/(M_m(j+1));
        disp([num2str(r_m(idx)), ' before and ', num2str(r_m_Dt(j+1)), ' after mixing'])
        
        r_m = r_m*0;                % reset r_m
        
    else
        r_m_Dt(j+1) = r_m_Dt(j);    %no impact, maintain same mantle Fe ratio
    end
end
        
if write == 1
    writematrix(t, fileOut, 'Sheet', sheetOut)
    writematrix(r_m_Dt, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)
end


figure(1);
hold on
box on
plot(t, r_m_Dt, "LineWidth", 1.5)
xlabel("Time (Myr)")
ylabel("Fe^{3+}/\SigmaFe Ratio")
hold off

