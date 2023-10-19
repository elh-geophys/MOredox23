% EL
% August 2022
% Updated 2023-09-29
%
% Produce Monte Carlo results with varying efficiency and depth for each
% impact during accretion.
% Uses getRainRatio.m

clear;

tic

% PARAMETERS TO CHANGE
model = 4;                  %accretion models: 4 = H04, 5 = N21  (using a continous model like 1-3 will take a LONG time :) )
r_0 = 0.004;                %initial Fe3+/sumFe
N = 1000;                   %number of MC samples
compSheet_early = 'EarthEarly';     %composition to use
compSheet_late = 'EarthLate';       %keep same as early if no change wanted
Tp_type = 'constant';               %chooose method to calculate Tp, either 'Pmo', 'U2Q', or 'constant'
    Tp_const = 3500;                %[K] Tp for 'constant' method
dP = 0.5e9;

sheetOut = 'H04_Tconst';            %sheet name to record data
fileOut = 'Rain_MC.xlsx';    % file name
write = 1;                   %1 or 0, to write to file

% ---------------------------------------------------------------------- %

% CONSTANTS
M_E_0 = 5.97e24;        %[kg] present day mass of Earth
M_c_0 = 1.88e24;        %[kg] present day mass of core
rho_imp = 5000;         %[kg/m^3] approximation based on weighted average (0.68Si + 0.32Fe)

% SET UP ACCRETION MODEL
[t, Accr_model] = getAccrModel(model);
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

% volume of core estimated by half the core impactor added
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

% determine the minimum pressure without impactor melt
P_min = zeros(1,length(M_imp));
for i = 1:length(M_imp)
    P_min(i) = get2LayerP(rho_m(i+1), rho_c_post(i), R_E_post(i), R_c_post(i), R_E_post_noimp(i));
end
P_min(isnan(P_min))=0;          %for some accretion models, where P=0 at t=0

[dMp_temp] = calcEqSi(M_c_imp, 'sph');

r_m_Dt = zeros(N,length(t));
switch Tp_type
    case 'Pmo'
        for i = 1:N
            if mod(i,50)==0
               disp(['Calculating trial ', num2str(i)]);
            end
            [r_m_Dt(i,:),~,~] = getRainRatio_Pmo(P_cmb, P_min, dP, r_0,Accr_model, ...
                dMp_temp, compSheet_early, compSheet_late, r_imp, rho_m, rho_c_post, M_m, M_m_imp, R_E_post, R_c_post);
        end
    case 'U2Q'
        for i = 1:N
            if mod(i,50)==0
               disp(['Calculating trial ', num2str(i)]);
            end
            [r_m_Dt(i,:),~,~] = getRainRatio_U2Q(P_cmb, P_min, dP, r_0, Accr_model, ...
                dMp_temp, compSheet_early, compSheet_late, r_imp, rho_m, rho_c_post, M_E, M_m, M_c, M_imp, M_m_imp, R_E, ...
                R_E_post, R_c_post, R_imp);
        end
    case 'constant'
        for i = 1:N
            if mod(i,50)==0
               disp(['Calculating trial ', num2str(i)]);
            end
            [r_m_Dt(i,:),~,~] = getRainRatio_constT(P_cmb, P_min, dP, r_0, Accr_model, ...
                dMp_temp, compSheet_early, compSheet_late, r_imp, rho_m, rho_c_post, M_m, M_m_imp, R_E_post, R_c_post, Tp_const);
        end
    otherwise
        disp('Could not calculate Tp based on input Tp_type')
end


% Determine the average and mid-50% and mid-90% 
r_prc = prctile(r_m_Dt,[0 1 5 25 50 75 95 99 100],1);
r_0p = r_prc(1,:);
r_1p = r_prc(2,:);
r_5p = r_prc(3,:);
r_25p = r_prc(4,:);
r_50p = r_prc(5,:);
r_75p = r_prc(6,:);
r_95p = r_prc(7,:);
r_99p = r_prc(8,:);
r_100p = r_prc(9,:);

r = [r_0p; r_1p; r_5p; r_25p; r_50p; r_75p; r_95p; r_99p; r_100p];

disp(r_50p(end))

if write == 1
    writematrix(t, fileOut, 'Sheet', sheetOut)
    writematrix(r, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)
end

toc

figure(2);
hold on
box on
plot(t, r)
xlabel("Time (Myr)")
ylabel("Fe^{3+}/\SigmaFe Ratio")
hold off

