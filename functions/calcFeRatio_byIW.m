% EL
% Sept 4, 2023
%
% Function to get the equilibrium (given dIW) Fe3+/sumFe value as a
% function of pressure, temperature (goetherm).
%
% INPUTS:   T           [K] geotherm temperature, scalar(?) or 1D array
%           P           [Pa] corresonding P for geotherm, scalar(?) or 1D array
%           dIW         [] dIW value for oxygen fugacity
%           PV_term     [] calculated int(PV)/RT from calcPV.m, scalar(?) or 1D array
%           OxiWts      [wt%] wt % for each compositional component
%           OxiMolW     [g/mol] mole weight for each component
%                 same order as compsheets: SiO2 TiO2 Al2O3 Cr2O3 FeO* MnO MgO NiO CaO Na2O K2O P2O5
%
% OUTPUTS:  FeRatio [] Fe3+/sumFe as a function of P, T, and dIW
%

function [FeRatio] = calcFeRatio_byIW(T, P, dIW, PV_term, OxiWts, OxiMolW)

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
    % See modPVcalc.m and calcPV.m to see how this is determined.
    PV_term = PV_term/log(10);

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

end

