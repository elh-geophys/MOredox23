% EL
% Feb 2023
% Updated 2023-09-06
%
% Eqn 21 of Hirschmann 2022
% Determine the Fe3+/sumFe ratio as a function of P,T for a given 
% metal-silicate equilibrium IW value.
% Use IW as a function of accretion time, based on Rubie+2011, Fig 3

clear;

%things to change: model, d_mo_factor (0 to 1), MoleWeights sheet, write
%files, t_idx to check a specific point
model = 4;          %accretion model to use, 4 = H04, 5 = N21
%rho_imp = rho_Si*0.68 + rho_Fe*0.32;            % weighted average for density of impactor
rho_imp = 5000;     %[kg/m^3] approximation based on weighted average (above)
M_E_0 = 5.97e24;    %[kg] present day mass of Earth
M_c_0 = 1.88e24;    %[kg] present day mass of core 

%note: only have this set up to do the Pmo method

d_mo_factor = 1;    %choose depth of melting, where 0 = surf and 1 = CMB

t_idx = 553;         %test a specific point in time

compSheet_early = 'EarthEarly';
compSheet_late = 'EarthLate';

write = 0;
writeFile = 'r_eq.xlsx';
writeSheet = 'N21';


% First 100 Myr accretion history of Earth
[t,Accr_model] = getAccrModel(model);
M_E = M_E_0 * Accr_model;
M_imp = M_E_0 * diff(Accr_model);

% note: during impact "n", earth mass is "n" pre-impact and "n+1" post-impact

% assume core takes up proportional mass of Earth
M_c = M_E * M_c_0/M_E_0;

% Rubie+2011, eqn 3
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

V_c = M_c./rho_c;                           %volume of core
V_m = (M_E - M_c)./rho_m;                   %volume of mantle
R_E = (3/(4*pi)*(V_m + V_c)).^(1/3);        %radius of Earth
R_c = (3/(4*pi)*(M_c./rho_c)).^(1/3);       %radius of core

%T&S Geodyanm Eqn 2.73, Pressure as a function radius from center with
P_cmb = zeros(1,length(R_c));
for i = 1:length(R_c)
    P_cmb(i) = get2LayerP(rho_m(i), rho_c(i), R_E(i), R_c(i), R_c(i));
end
P_cmb(isnan(P_cmb))=0;

% choose impactor MO depth
P_mo = P_cmb*d_mo_factor;
P_max = round(P_mo,-9);

Tp = getMOTp(P_mo/1e9);
Tad = zeros((P_max(end)/1e9*2)+1,length(P_max));      %set up Tad

%determine adiabat, from code by J.Korenaga
%this will get adiabats for whole accretion time
P = ((0:1e9:(P_max(1)*2))/2)';           %to go by 0.5 GPa
if size(P)==1
    P = [0, 0.5e9];
end

Tad_temp = getMOAdiabat(Tp(1),P);
Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
Tad_temp2 = cat(1, Tad_temp, Tad_empty);
Tad = Tad_temp2;

for i = 2:length(P_max)
    if M_imp(i-1) > 0
        P = ((0:1e9:(P_max(i)*2))/2)';           %to go by 0.5 GPa
    
        Tad_temp = getMOAdiabat(Tp(i),P);
        Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);
        
    else
        Tad_temp = Tad(:,i-1);
        Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);
    end
end

T = Tad(:,t_idx);           %choose one point in time
T(isnan(T))=[];

P = P(1:length(T));

% read desired geotherm
% geodata = readmatrix('\db\geotherms.xlsx', 'Sheet', geoSheet);
% T = geodata(:,4);
% P = geodata(:,2);


% dIW value considered metal-silicate equilibrium, see Rubie+2011
F = Accr_model(t_idx);
dIW = 3*F-5;

% Fitting parameters from Table 1 of Hirschmann 2021
R = 8.314;
T0 = 1673.15;
a = 0.1917;
b = -1.961;
c = 4158.1;
Cp = 33.25;
y1 = -520.46;   %SiO2    param from Borisov+2018, Eqn 20 of H22
y2 = -185.37;   %TiO2
y3 = 494.39;    %MgO
y4 = 1838.34;   %CaO
y5 = 2888.48;   %Na2O
y6 = 3473.68;   %K2O
y7 = -4473.6;   %PO5
y8 = -1245.09;  %SiO2*Al2O3
y9 = -1156.86;  %SiO2*MgO

% read %wt and mol wt data
% SiO2 TiO2 Al2O3 Cr2O3 FeO* MnO MgO NiO CaO Na2O K2O P2O5
% data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet);
% OxiWts = data(:,2);
% OxiMolW = data(:,3);

data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet_early);
OxiWts_early = data(:,2);
OxiMolW = data(:,3);

if compSheet_late == 0
    compSheet_late = compSheet_early;       %no change
    data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet_late);
    OxiWts_late = data(:,2);
else
    data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet_late);
    OxiWts_late = data(:,2);
end

slope = OxiWts_late - OxiWts_early;
OxiWts = slope.*Accr_model + OxiWts_early;

tempMol = OxiWts./OxiMolW;
MolSum = sum(tempMol);
Mol = tempMol./MolSum;

Mol = Mol(:,t_idx);

% Activity ratio term
ActRatios_term = (1./T) .* (y1*Mol(1) + y2*Mol(2) + y3*Mol(7) + y4*Mol(9) +...
    y5*Mol(10) + y6*Mol(11) + y7*Mol(12) + y8*Mol(1)*Mol(3) + y9*Mol(1)*Mol(7));

% Gibbs energy term
Gb_term = Cp/(R*log(10)) * (1 - T0./T - log(T./T0));

% Work term
% Using pre-calculated PV term for the geotherm
% This is = integration/(R*T) so still need to 1/ln(10) for H22 Eqn 21
% See modPVcalc.m to see how this is determined.
% PVtemp1 = readmatrix('\db\PVcalc.xlsx');
% PVtemp2 = PVtemp1(:,PVcol);
% PV_term = PVtemp2/log(10);
PV_term_temp = calcPV(T,P)/log(10);
PV_empty = NaN(length(P)-length(PV_term_temp),1);
PV_term = cat(1,PV_term_temp,PV_empty);

% Fugacity term
IW = zeros(length(P),1);
for i = 1:length(P)
    if P(i) <= 100e9
        IW(i) = getIW_H21(P(i),T(i));
        i_last = i;
    else     %linear extrap past 100GPa
        IW(i) = interp1(P(1:i_last), IW(1:i_last), P(i), 'linear', 'extrap');
    end
end

logfO2 = dIW + IW;
fO2_term = a*logfO2 + b + c./T;

% Eqn 21 of Hirschmann 2022, as log10(FeO1.5/FeO)
logFeOxiRatio = fO2_term - Gb_term - PV_term + ActRatios_term;

%convert to Fe3/Fe_total
FeRatio = 1./(1+(1./10.^(logFeOxiRatio)));


FeRatio_test = calcFeRatio(T, P, F, PV_term*log(10), compSheet_early, compSheet_late);

if write == 1
    writematrix(P, writeFile, 'Sheet', writeSheet, 'Range', 'A2')
    writematrix(t, writeFile, 'Sheet', writeSheet, 'Range', 'B1')
    writematrix(FeRatio, writeFile, 'Sheet', writeSheet, 'Range', 'B2') 
end

figure(1);
hold on
box on
x = [0.03 0.045];       %shaded region
y = [135 135];
area(x,y, 'FaceColor', [0.75 0.75 0.75])
plot(FeRatio, P/1e9, 'b-', "LineWidth", 1.5);
plot(FeRatio_test, P/1e9, 'r--', "LineWidth", 1.5);
xlabel('Fe^{3+}/\SigmaFe Ratio')
ylabel('Pressure (GPa)')
ax = gca;
ax.YDir = 'reverse';

figure(2);
plot(T, P/1e9)
xlabel("Temperature (K)")
ylabel("Pressure (GPa)")

figure(3);
hold on
plot(fO2_term,P/1e9)
plot(Gb_term,P/1e9)
plot(PV_term,P/1e9)
plot(ActRatios_term,P/1e9)
plot(logFeOxiRatio,P/1e9, "LineWidth",2)
xlim([-3 4])
legend("fO2", "Gb", "PV", "ActRatio", "logFe3/Fe2", "Location", "southeast")
