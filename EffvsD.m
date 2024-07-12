% EL
% Original: August 2022
% Updated: 07-02-2024
%
% Determine Fe3/sumFe in MO as a function of given depth (P_GI) and given
% efficiency (eff) for the last giant impact (GI)
%
% Creates data for Figure 3 in MO redox manuscript. Use this data in the
% corresponding ***_fig.m file to output the figure. This script outputs a
% figure, but it is unpolished.

clear;

% PARAMETERS
model = 5;                          %accretion models, 4 = H04, 5 = N21
r_0_idx = 2;                        %index for r_0, 1=0th, 2=1st, 3=5th... 8=99th, 9=100th
compSheet_earth = 'N21_E';     %sheet in Compositions.xlsx to use for composition
compSheet_imp = 'N21_imp';
Tp_type = 'Pmo';               %constant or Pmo
dP = 0.5e9;                         %[Pa] increments of P for layer to do metal rain calculation

sheet_nomix = 'N21_1st_nomix';     %sheet name to record data
sheet_mix = 'N21_1st_mix';
dataSheet = 'N21_data';
fileOut = 'Rain_EffvsD_Pmo.xlsx';       % file name
write = 1;                          %1 to write, else to not write to file

%CHOOSE YOUR Fe3/sumFe VALUE BEFORE GI
%                 [0th    1st    5th    25th   50th   75th   95th   99th   100th]
if model == 4   % column UF in Rain_MC.xlsx
    %H04
    switch Tp_type
        case 'constant'
            r_0 = [0.0695,0.0781,0.0878,0.1009,0.1101,0.1189,0.1271,0.1299,0.1314];    %Tconst
        case 'Pmo'
            r_0 = [0.0365,0.0523,0.0709,0.0963,0.1156,0.1316,0.1470,0.1510,0.1538];    %Pmo
    end
elseif model == 5    % column OT
    %N21
    switch Tp_type
        case 'constant'
            r_0 = [0.0569,0.0644,0.0804,0.0953,0.1047,0.1108,0.1210,0.1267,0.1289];    %Tconst
        case 'Pmo'
            r_0 = [0.0311,0.0367,0.0426,0.0615,0.0793,0.1052,0.1285,0.1338,0.1362];    %Pmo
    end
else
    disp('No r_0 values for model chosen')
end

% ---------------------------------------------------------------------- %

% READ DATA SHEETS
PV_data = readmatrix('/db/PVcalc.xlsx');
Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');
CompEarth_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_earth);
CompImp_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_imp);
MolW_byM_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Emantle', 'Range', 'B2:B13');
MolW_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Emantle', 'Range', 'C2:C13');

% CONSTANTS
M_E_0 = 5.97e24;        %[kg] present day mass of Earth
rho_imp = 5000;         %[kg/m^3] approximation based on weighted average (0.68Si + 0.32Fe)

% SET UP ACCRETION MODEL
[t, Accr_model] = getAccrModel(model);
GI_idx = find(diff(Accr_model),1,'last');         %we only want for GI

M_E = M_E_0 * Accr_model(GI_idx);
M_imp = M_E_0 * (Accr_model(GI_idx+1)-Accr_model(GI_idx));

FeO_E = CompEarth_data(7,GI_idx);

% core/mantle evolution
M_c_ratio = CompEarth_data(2,:);
M_c = M_E * M_c_ratio(GI_idx);    %core
M_m = M_E - M_c;                  %mantle
M_m_post = M_E_0*Accr_model(GI_idx+1) - M_E*M_c_ratio(GI_idx+1);

% Rubie+2011, approximations for mantle and core densities
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

M_c_imp = M_E_0*Accr_model(GI_idx+1)*M_c_ratio(GI_idx+1) - M_c;
M_m_imp = M_imp - M_c_imp;

%estimated Fe3+/sumFe for impactor, based on Rubie+2011, Supp. Table 3a
r_imp = CompImp_data(5,GI_idx);
FeO_imp = CompImp_data(4,GI_idx);

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
            otherwise
                disp('Could not calculate Tp based on input Tp_type')
                Tp = NaN;
                break               % ends the if statement
        end
        
        Tad = getMOAdiabat(Tp, P, Adiabat_data);
        
        % DETERMINE PV INTEGRATION AS INT(PV)/RT
        PV = calcPV(Tad,P,PV_data);
        
        % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
        [r_eq,~] = calcFeRatio(Tad, P, PV, CompEarth_data(3:end,GI_idx+1), MolW_data, MolW_byM_data);
     
        % METAL RAIN CALCULATION FOR THIS TIME STEP
        R_mo = interp1(P_check, R, P_GI(j));
        M_mo = rho_m(GI_idx+1) * 4/3*pi*(R_E_post^3 - R_mo^3);        %approximate mass of spherical shell MO
    
        % initial Fe ratio of mantle after mixing previous silicate with impactor silicate
        r_m(1) = ((M_mo-M_m_imp)*FeO_E*r_0(r_0_idx)+M_m_imp*FeO_imp*r_imp)/((M_mo-M_m_imp)*FeO_E + M_m_imp*FeO_imp);
        FeO_mo = ((M_mo-M_m_imp)*FeO_E + M_m_imp*FeO_imp)/M_mo;
        
        dMp = min(dMp_temp*eff(k), M_mo);                               %silicate mass eq'd with given efficiency
        
        for i = 1:length(P)        % droplets falling through mantle time
            % ratio changes as the droplets fall through mantle layers until base of MO
            % note r_m(i+1) is at P(i) location
            r_m(i+1) = (dMp*r_eq(i) + (M_mo-dMp)*r_m(i))/M_mo;
            idx = i+1;
        end
        
        % final r from mixing redox'd MO with whole mantle
        r_m_mix = (M_mo*FeO_mo*r_m(idx)+(M_m_post-M_mo)*FeO_E*r_0(r_0_idx))/(M_mo*FeO_mo + (M_m_post-M_mo)*FeO_E);
        
        FeO_f = (M_mo*FeO_mo+(M_m_post-M_mo)*FeO_E)/(M_m_post);

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


%contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16];
c0 = [0.02 0.06];

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

map = [0.82 0.82 0.82; 1 1 1];   
%map = [1 1 1];                                      
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


