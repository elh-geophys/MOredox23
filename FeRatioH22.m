% EL
% Feb 2023
% Updated 2023-03-01
%
% Eqn 21 of Hirschmann 2022
% Determine the Fe3+/sumFe ratio as a function of P,T for a given 
% metal-silicate equilibrium IW value.

clear;

%things to change: geotherm data sheet, PVcol, MoleWeights sheet, write column
geoSheet = '4500K';
PVcol = 7;              %2=2000K, 3=2500K, 4=3000K, 5=3500K, etc.
compSheet = 'Armstrong19';
writeSheet = 'Armstrong19';
writeCol = 'G1';        %B1=2000K, C1=2500K, D1=3000K, etc.

% dIW value considered metal-silicate equilibrium
dIW = -2;

% read desired geotherm
geodata = readmatrix('\db\geotherms.xlsx', 'Sheet', geoSheet);
T = geodata(:,4);
P = geodata(:,2);

% Fitting parameters from Table 1
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
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet);
OxiWts = data(:,2);
OxiMolW = data(:,3);

tempMol = OxiWts./OxiMolW;
MolSum = sum(tempMol);
Mol = tempMol/MolSum;

% Activity ratio term
ActRatios_term = (1./T) .* (y1*Mol(1) + y2*Mol(2) + y3*Mol(7) + y4*Mol(9) +...
    y5*Mol(10) + y6*Mol(11) + y7*Mol(12) + y8*Mol(1)*Mol(3) + y9*Mol(1)*Mol(7));

% Gibbs energy term
Gb_term = Cp/(R*log(10)) * (1 - T0./T - log(T./T0));

% Work term
% Using pre-calculated PV term for the geotherm
% This is = integration/(R*T) so still need to 1/ln(10) for H22 Eqn 21
% See modPVcalc.m to see how this is determined.
PVtemp1 = readmatrix('\db\PVcalc.xlsx');
PVtemp2 = PVtemp1(:,PVcol);
PV_term = PVtemp2/log(10);

% Fugacity term
IW = zeros(length(P),1);
for i = 1:length(P)
    if P(i) <= 100
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

writematrix(P, 'r_eq.xlsx', 'Sheet', writeSheet, 'Range', 'A1')
writematrix(FeRatio, 'r_eq.xlsx', 'Sheet', writeSheet, 'Range', writeCol) 

figure(1);
hold on
box on
yyaxis left
colororder('default')
x = [0.03 0.045];       %shaded region
y = [135 135];
area(x,y, 'FaceColor', [0.75 0.75 0.75])
r5 = plot(FeRatio, P, "Color", "#77AC30", "LineWidth", 1.5);
xlabel('Fe^{3+}/\SigmaFe Ratio')
ylabel('Pressure (GPa)')
