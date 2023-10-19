% EL
% August 2022
% Updated 2023-09-29
%
% Determine Fe3/Fe in MO of given depth (P_mo) and given efficiency (eff)
% for late accretion (1%)


clear;

% PARAMETERS
compSheet_early = 'EarthEarly';     %sheet in MoleWeights.xlsx to use for composition
compSheet_late = 'EarthLate';
r_imp = 0.0122;                     %Fe3+/sumFe of impactor silicate
Tp_type = 'Pmo';                    %constant, Pmo, or U2Q
    T0 = 1613;                          %for U2Q method
dP = 0.5e9;    

%Fe3/sumFe value AFTER GI
%Tconst method
    %r_0 = 0.0862;      %1st H04
    %r_0 = 0.0949;      %5th H04  
    %r_0 = 0.1078;      %25th H04
    %r_0 = 0.1177;      %median from H04 modeling
    %r_0 = 0.1134;      %median from N21 modeling
    %r_0 = 0.1030;      %25th N21
    %r_0 = 0.0949;      %5th N21
%Pmo method
     r_0 = 0.0718;       %1st H04  
%     r_0 = 0.0845;       %5th H04
%     r_0 = 0.1085;       %25th H04
%     r_0 = 0.1272;       %median H04
%     r_0 = 0.0918;       %5th N21
%     r_0 = 0.1040;       %25th N21
%     r_0 = 0.1194;       %median N21
%U2Q method
%     r_0 = 0.0507;       %5th H04
%     r_0 = 0.0641;       %25th H04
%     r_0 = 0.0808;       %median H04
%     r_0 = 0.0901;       %5th N21
%     r_0 = 0.1033;       %25th N21
%     r_0 = 0.1189;       %median N21

sheet_nomix = 'H04_1st_nomix';     %sheet name to record data
sheet_mix = 'H04_1st_mix';
dataSheet = 'data';
fileOut = 'Rain_EffvsD_late_Pmo.xlsx';               % file name
write = 0;

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
Accr_model = [0.99, 1];

M_E = Accr_model(1) * M_E_0;         % mass post-GI
M_imp = diff(Accr_model) * M_E_0;       % mass accreted in late accretion

% assume core and mantle take up proportional mass of Earth
M_c = M_E * M_c_0/M_E_0;    %core
M_m = M_E - M_c;            %mantle
M_m_post = M_E_0 - M_c_0;

% Rubie+2011, approximations for mantle and core densities
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

V_c = M_c/rho_c(1);                           %volume of core
V_m = M_m/rho_m(1);                           %volume of mantle
R_E = (3/(4*pi)*(V_m + V_c))^(1/3);           %radius of Earth
R_imp = (3/(4*pi)*M_imp/rho_imp)^(1/3);       %radius of impactor 

Fe_ratio = 0.321;               % [] weight % of iron on Earth/impactor
Si_ratio = 1-Fe_ratio;          % [] weight % of silicate on Earth/impactor
M_c_imp = M_imp*Fe_ratio;       % [kg] approximate proportion metal mass of impactor
M_m_imp = M_imp*Si_ratio;       % [kg] approximate proportion silicate mass of impactor

% Use half impactor core mass between pre- and post-impact to determine R_c at impact
% Use entire mantle mass for chemical mixing post-impact
M_c_post = M_c+M_c_imp/2;                         %core mass with half of the impactor 
rho_c_post = (rho_c(1)+rho_c(2))/2;                  %average core density between pre and post impact
R_c_post = (3/(4*pi)*(M_c_post/rho_c_post))^(1/3);   %radius of core post impact, Earth + half of impactor

% radius of Earth w/ all impactor Si (M_m post-impact) and 1/2 impactor Fe
R_E_post = (3/(4*pi)*(M_m_post/rho_m(2) + M_c_post/rho_c_post))^(1/3);

% radius of Earth w/ no impactor Si (M_m pre-impact), but 1/2 impactor Fe;
% upper bound for MO radius with all impactor as melt
R_E_post_noimp = (3/(4*pi)*(M_m/rho_m(2) + M_c_post/rho_c_post))^(1/3);

%T&S Geodyanm Eqn 2.73, Pressure as a function radius from center
% determine pressure at CMB
R = linspace(R_c_post,R_E_post,1000);
P_check = get2LayerP(rho_m(2), rho_c_post, R_E_post, R_c_post, R);
P_cmb = P_check(1);
P_max = P_cmb-mod(P_cmb,dP);

% determine the minimum pressure without impactor melt
P_min = get2LayerP(rho_m(2), rho_c_post, R_E_post, R_c_post, R_E_post_noimp);
P_min = P_min+dP-mod(P_min,dP);

[dMp_temp] = calcEqSi(M_c_imp, 'sph');

% DETERMINE POST-IMPACT FE RATIO
eff = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]/100;
P_late = P_min:dP:P_max;

r_m = zeros(length(P_late),1);
effvsd_nomix = zeros(length(P_late), length(eff));
effvsd_mix = zeros(length(P_late), length(eff));

for k = 1:length(eff)          % for each efficiency
    disp(['Calculating for e = ', num2str(eff(k))])
    
    for j = 1:length(P_late)     % for each depth
    
        P = (0:dP:P_late(j))';            % initiate a P array for this GI depth

        % GET MO ADIABAT
        switch Tp_type
            case 'Pmo'
                Tp = getMOTp(P_late(j)/1e9, Adiabat_data);
            case 'constant'
                Tp = 3500;
            case 'U2Q'
                Tp_mo_base = getMOTp(P_late(j)/1e9, Adiabat_data);
                
                Tf_test = calcTforU2Q(T0,Cp,M_E,M_c,M_imp,R_E,R_imp,Si_ratio,0.2);
                if Tf_test > Tp_mo_base
                    Tp = Tp_mo_base;
                else
                    Tp = Tf_test;
                end
            otherwise
                disp('Could not calculate Tp based on input Tp_type')
                Tp = NaN;
                break               % ends the if statement
        end
        
        Tad = getMOAdiabat(Tp, P, Adiabat_data);
        
        % DETERMINE PV INTEGRATION AS INT(PV)/RT
        PV = calcPV(Tad,P, PV_data);
        
        % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
        [r_eq,~] = calcFeRatio(Tad, P, Accr_model(2), PV, CompEarly_data, CompLate_data);
     
        % METAL RAIN CALCULATION FOR THIS TIME STEP
        R_mo = interp1(P_check, R, P_late(j));
        M_mo = rho_m(2) * 4/3*pi*(R_E_post^3 - R_mo^3);        %approximate mass of spherical shell MO
    
        % initial Fe ratio of mantle after mixing previous silicate with impactor silicate
        r_m(1) = ((M_mo-M_m_imp)*r_0+M_m_imp*r_imp)/(M_mo);
        
        dMp = min(dMp_temp*eff(k), M_mo);    %silicate mass eq'd with given efficiency
        
        for i = 1:length(P)        % droplets falling through mantle time
            % ratio changes as the droplets fall through mantle layers until base of MO
            % note r_m(i+1) is at P(i) location
            r_m(i+1) = (dMp*r_eq(i) + (M_mo-dMp)*r_m(i))/(M_mo);
            idx = i+1;
        end
        
        % final r from mixing redox'd MO with whole mantle
        r_m_mix = (M_mo*r_m(idx) + (M_m_post-M_mo)*r_0)/M_m_post;

        effvsd_nomix(j,k) = r_m(idx);
        effvsd_mix(j,k) = r_m_mix;
        
        r_m = r_m*0;                % reset r_m
    end

end

if write == 1
    writematrix(effvsd_nomix, fileOut, 'Sheet', sheet_nomix)
    writematrix(effvsd_mix, fileOut, 'Sheet', sheet_mix)
    writematrix(eff, fileOut, 'Sheet', dataSheet)
    writematrix(P_late, fileOut, 'Sheet', dataSheet, 'Range', 'A2')
end


% THIS MAY OR MAY NOT WORK :)

c = linspace(0.005, 0.1, 20);
ct = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10];
c0 = [0.02 0.03];

figure('Position', [200 200 1000 400]);
subplot(1,2,1)
hold on
box on
[C,h] = contour(P_late/1e9, log10(eff), effvsd_mix'-0.35/8.1, c, 'LineWidth', 2);
clabel(C,h,ct, 'LabelSpacing', 200)
xlabel('Pressure (GPa)')
ylabel('log(Efficiency)')
xlim([P_min/1e9 100])
title("Whole Mantle (with post-mixing)")

subplot(1,2,2)
hold on
box on
colormap cool
colormap(flipud(colormap));
contour(P_late/1e9, log10(eff), effvsd_nomix'-0.35/8.1, c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 300)
xlabel('Pressure (GPa)')
ylabel('log(Efficiency)')
xlim([P_min/1e9 100])
ylim([-3 0])
title("Top MO only (no post-mixing)")
