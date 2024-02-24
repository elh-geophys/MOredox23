% EL
% Original: August 2022
% Updated: 02-22-2024
%
% Determine Fe3/sumFe in MO as a function of given depth (P_GI) and given
% efficiency (eff) for the last giant impact (GI)

clear;

% PARAMETERS
model = 5;                          %accretion models, 4 = H04, 5 = N21
r_0_idx = 4;                        %index for r_0, 1=0th, 2=1st, 3=5th... 8=99th, 9=100th
compSheet_earth = 'N21_E';     %sheet in Compositions.xlsx to use for composition
compSheet_imp = 'N21_imp';
Tp_type = 'Pmo';               %constant, Pmo, or U2Q method
    T0 = 1613;                          %initial temp for U2Q
dP = 0.5e9;                         %[Pa] increments of P for layer to do metal rain calculation

sheet_nomix = 'N21_25th_nomix';     %sheet name to record data
sheet_mix = 'N21_25th_mix';
dataSheet = 'N21_data';
fileOut = 'Rain_EffvsD_Pmo.xlsx';       % file name
write = 1;                          %1 to write, else to not write to file

% ---------------------------------------------------------------------- %

% READ DATA SHEETS
PV_data = readmatrix('/db/PVcalc.xlsx');
Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');
CompEarth_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_earth);
CompImp_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_imp);
MolW_byM_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11', 'Range', 'B2:B13');
MolW_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11', 'Range', 'C2:C13');

%CHOOSE YOUR Fe3/sumFe VALUE BEFORE GI
%         [0th    1st    5th    25th   50th   75th   95th   99th   100th]
if model == 4   % column UG in Rain_MC.xlsx
    %H04
    %r_0 = [0.0651,0.0794,0.0867,0.1001,0.1085,0.1170,0.1263,0.1297,0.1314];     %Tconst
    r_0 = [0.0415,0.0533,0.0680,0.0950,0.1146,0.1319,0.1474,0.1532,0.1555];     %Pmo
elseif model == 5    % column OU
    %N21
    %r_0 = [0.0317,0.0381,0.0448,0.0533,0.0591,0.0629,0.0690,0.0724,0.0742];     %Tconst
    r_0 = [0.0160,0.0193,0.0230,0.0330,0.0440,0.0583,0.0728,0.0765,0.0785];     %Pmo
else
    disp('No r_0 values for model chosen')
end

% CONSTANTS
Cp = 1e3;               %[J/kgK] specific heat 
M_E_0 = 5.97e24;        %[kg] present day mass of Earth
M_c_0 = 1.88e24;        %[kg] present day mass of core
rho_imp = 5000;         %[kg/m^3] approximation based on weighted average (0.68Si + 0.32Fe)

% SET UP ACCRETION MODEL
[t, Accr_model] = getAccrModel(model);
GI_idx = find(diff(Accr_model),1,'last');         %we only want for GI

M_E = M_E_0 * Accr_model(GI_idx);
M_imp = M_E_0 * (Accr_model(GI_idx+1)-Accr_model(GI_idx));

FeO_E = CompEarth_data(5,GI_idx);

% assume core and mantle take up proportional mass of Earth
M_c = M_E * M_c_0/M_E_0;    %core
M_m = M_E - M_c;            %mantle
M_m_post = M_E_0*Accr_model(GI_idx+1) - M_c_0*Accr_model(GI_idx+1);

% Rubie+2011, approximations for mantle and core densities
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

V_c = M_c/rho_c(GI_idx);                           %volume of core
V_m = M_m/rho_m(GI_idx);                           %volume of mantle
R_E = (3/(4*pi)*(V_m + V_c))^(1/3);                %radius of Earth
R_imp = (3/(4*pi)*M_imp/rho_imp)^(1/3);            %radius of impactor 

Fe_ratio = 0.321;               % [] weight % of iron on Earth/impactor
Si_ratio = 1-Fe_ratio;          % [] weight % of silicate on Earth/impactor
M_c_imp = M_imp*Fe_ratio;       % [kg] approximate proportion metal mass of impactor
M_m_imp = M_imp*Si_ratio;       % [kg] approximate proportion silicate mass of impactor

%estimated Fe3+/sumFe for impactor, based on Rubie+2011, Supp. Table 3a
r_imp = CompImp_data(3,GI_idx)+0.35/8.2;
FeO_imp = CompImp_data(2,GI_idx);

% Use half impactor core mass between pre- and post-impact to determine R_c at impact
% Use entire mantle mass for chemical mixing post-impact
M_c_post = M_c+M_c_imp/2;                            %core mass with half of the impactor 
rho_c_post = (rho_c(GI_idx)+rho_c(GI_idx+1))/2;      %average core density between pre and post impact
R_c_post = (3/(4*pi)*(M_c_post/rho_c_post))^(1/3);   %radius of core post impact, Earth + half of impactor

% radius of Earth w/ all impactor Si (M_m post-impact) and 1/2 impactor Fe
R_E_post = (3/(4*pi)*(M_m_post/rho_m(GI_idx+1) + M_c_post/rho_c_post))^(1/3);

% radius of Earth w/ no impactor Si (M_m pre-impact), but 1/2 impactor Fe;
% upper bound for MO radius with all impactor as melt
R_E_post_noimp = (3/(4*pi)*(M_m/rho_m(GI_idx+1) + M_c_post/rho_c_post))^(1/3);

%T&S Geodyanm Eqn 2.73, Pressure as a function radius from center
% determine pressure at CMB
R = linspace(R_c_post,R_E_post,1000);
P_check = get2LayerP(rho_m(GI_idx+1), rho_c_post, R_E_post, R_c_post, R);
P_cmb = P_check(1);
P_max = P_cmb-mod(P_cmb,dP);

% determine the minimum pressure without impactor melt
P_min = get2LayerP(rho_m(GI_idx+1), rho_c_post, R_E_post, R_c_post, R_E_post_noimp);
P_min = P_min+dP-mod(P_min,dP);

[dMp_temp] = calcEqSi(M_c_imp, 'sph');

% DETERMINE POST-GI FE RATIO
eff = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]/100;
P_GI = P_min:dP:P_max;

r_m = zeros(length(P_GI)+1,1);
effvsd_nomix = zeros(length(P_GI), length(eff));
effvsd_mix = zeros(length(P_GI), length(eff));

for k = 1:length(eff)          % for each efficiency
    disp(['Calculating for e = ', num2str(eff(k))])
    
    for j = 1:length(P_GI)     % for each depth

        P = (0:dP:P_GI(j))';            % initiate a P array for this GI depth
        
        % GET MO ADIABAT
        switch Tp_type
            case 'Pmo'
                Tp = getMOTp(P_GI(j)/1e9, Adiabat_data);
            case 'constant'
                Tp = 3500;
            case 'U2Q'
                Tp_mo_base = getMOTp(P_GI(j)/1e9, Adiabat_data);
                
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
        PV = calcPV(Tad,P,PV_data);
        
        % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
        [r_eq,~] = calcFeRatio(Tad, P, PV, CompEarth_data(:,GI_idx+1), MolW_data, MolW_byM_data);
     
        % METAL RAIN CALCULATION FOR THIS TIME STEP
        R_mo = interp1(P_check, R, P_GI(j));
        M_mo = rho_m(GI_idx+1) * 4/3*pi*(R_E_post^3 - R_mo^3);        %approximate mass of spherical shell MO
    
        % initial Fe ratio of mantle after mixing previous silicate with impactor silicate
        %r_m(1) = ((M_mo-M_m_imp)*r_0(r_0_idx)+M_m_imp*r_imp)/(M_mo);
        r_m(1) = ((M_mo-M_m_imp)*FeO_E*r_0(r_0_idx)+M_m_imp*FeO_imp*r_imp)/((M_mo-M_m_imp)*FeO_E + M_m_imp*FeO_imp);
        FeO_mo = ((M_mo-M_m_imp)*FeO_E + M_m_imp*FeO_imp)/M_mo;
        %disp(['Initial Fe3+/sumFe MO = ', num2str(r_m(1)), ' and FeO of MO = ', num2str(FeO_mo)]);
        
        dMp = min(dMp_temp*eff(k), M_mo);                               %silicate mass eq'd with given efficiency
        
        for i = 1:length(P)        % droplets falling through mantle time
            % ratio changes as the droplets fall through mantle layers until base of MO
            % note r_m(i+1) is at P(i) location
            r_m(i+1) = (dMp*r_eq(i) + (M_mo-dMp)*r_m(i))/M_mo;
            idx = i+1;
        end
        
        % final r from mixing redox'd MO with whole mantle
        %r_m_mix = (M_mo*r_m(idx) + (M_m_post-M_mo)*r_0(r_0_idx))/M_m_post;
        r_m_mix = (M_mo*FeO_mo*r_m(idx)+(M_m-M_mo)*FeO_E*r_0(r_0_idx))/(M_mo*FeO_mo + (M_m_post-M_mo)*FeO_E);
        %disp([num2str(r_m(idx)), ' before and ', num2str(r_m_mix), ' after mixing'])

        FeO_f = (M_mo*FeO_mo+(M_m_post-M_mo)*FeO_E)/(M_m_post);
        %disp(['FeO wt% = ', num2str(FeO_f)])

        effvsd_nomix(j,k) = r_m(idx);
        effvsd_mix(j,k) = r_m_mix;
        
        r_m = r_m*0;                % reset r_m
        
    end

end

if write == 1
    writematrix(effvsd_nomix, fileOut, 'Sheet', sheet_nomix)
    writematrix(effvsd_mix, fileOut, 'Sheet', sheet_mix)
    writematrix(eff, fileOut, 'Sheet', dataSheet)
    writematrix(P_GI, fileOut, 'Sheet', dataSheet, 'WriteMode', 'append')
end

%range for post-Cr oxidation = modern day mantle FeO*
r_low_f = 0.02;
r_high_f = 0.06;

%range for pre-Cr oxidation with 8.1% FeO* from Deng20 composition
r_low_0 = r_low_f + 0.35/8.2;
r_high_0 = r_high_f + 0.35/8.2;

%contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16];
c0 = [0.02 0.03 0.04 0.05 0.06];

% THIS MAY OR MAY NOT WORK :)
figure('Position', [200 200 500 400]);
ax1 = axes;
contourf(ax1, P_GI/1e9, log10(eff), (effvsd_mix'-0.35/8.2), c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, P_GI/1e9, log10(eff), (effvsd_nomix'-0.35/8.2), c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

map = [0.88 0.88 0.88; 0.80 0.80 0.80; 0.72 0.72 0.72; 1 1 1; 1 1 1];   %use depending on # of contours
%map = [0.72 0.72 0.72; 1 1 1; 1 1 1];                                      
colormap(ax1, map)
colormap(ax2, flipud(colormap(ax2,cool)));

ax1.Box = 'on';
ax1.XLabel.String = 'P (GPa)';
ax1.YLabel.String = 'Efficiency';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'0.1%' '' '' '' '' '' '' '' '' ...
    '1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};


