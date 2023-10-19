% EL
% Feb 2023
%
% This function gets the log10(fO2) value for a given P, T, and Fe3+/sumFe
% Based on Hirschmann2022 Eqn 21
%
% Inputs:  P, T, Fe3+/sumFe ("r"), sheet for composition (in
% MoleWeights.xlsx)
% Output: log10(fO2) compared to IW (DeltaIW), and logfO2

function [logfO2vsIW, logfO2] = getfO2_H22(P, T, PV_term, r, compSheet)

    % read %wt and mol wt data
    % SiO2 TiO2 Al2O3 Cr2O3 FeO* MnO MgO NiO CaO Na2O K2O P2O5
    data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', compSheet);
    OxiWts = data(:,2);
    OxiMolW = data(:,3);

    tempMol = OxiWts./OxiMolW;
    MolSum = sum(tempMol);
    Mol = tempMol/MolSum;

    % Fitting parameters from H22, Table 1
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

    % Activity ratio term
    ActRatios_term = (1./T) .* (y1*Mol(1) + y2*Mol(2) + y3*Mol(7) + y4*Mol(9) +...
        y5*Mol(10) + y6*Mol(11) + y7*Mol(12) + y8*Mol(1)*Mol(3) + y9*Mol(1)*Mol(7));

    % Gibbs energy term
    Gb_term = Cp/(R*log(10)) * (1 - T0./T - log(T./T0));

    % Work term
    % Using pre-calculated PV term for the geotherm
    % This is = integration/(R*T) so still need to 1/ln(10) for H22 Eqn 21
    % see modPVcalc.m to see how this is calculated from Deng2020 EOS
%     Tinit = T(1);
%     name = "PV_"+Tinit+"KTp";
%     table = readtable('\db\PVcalc_old.xlsx','ReadVariableNames',true);
%     PVtemp = table{:,name};
%     PV_term = PVtemp/log(10);

    % FeRatio term
    % Note that in H22 Eqn 21, the FeRatio is given as log10(Fe3+/Fe2+), so
    % have to calculate from r_eq as Fe3+/sumFe
    FeOxiRatio = r ./ (1-r);
    logFeOxiRatio = log10(FeOxiRatio);

    % Fugacity term
    fO2_term = logFeOxiRatio + Gb_term + PV_term - ActRatios_term;
    logfO2 = (1/a)*(fO2_term - b - c./T);

    % compare to IW
    IW = zeros(length(P),1);
    for i = 1:length(P)
        if P(i) <= 100e9
            IW(i) = getIW_H21(P(i),T(i));
            i_last = i;
        else     %linear extrap past 100GPa
            IW(i) = interp1(P(1:i_last), IW(1:i_last), P(i), 'linear', 'extrap');
        end
    end

    logfO2vsIW = logfO2 - IW;

end 
