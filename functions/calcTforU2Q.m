% EL
% Aug 31, 2023
%
% Function to calculate the GPE (U) between Earth and impactor and the
% corresponding temperature increase from heating (Q).
%
% INPUTS:   T0       [K] initial temperature of mantle, e.g. 1613K
%           Cp       [J/kgK] heat capacity of mantle/MO
%           M_E      [kg] mass of Earth
%           M_c      [kg] mass of Earth's core
%           M_imp    [kg] mass of impactor
%           R_E      [m] radius of Earth
%           R_imp    [m] radius of impactor
%           Si_ratio [] fraction of impactor that is silicate mantle
%           h        [] factor for U that goes to Q
%
% OUTPUTS:  Tf      [K] final temperature of silicate mantle
%

function [Tf] = calcTforU2Q(T0, Cp, M_E, M_c, M_imp, R_E, R_imp, Si_ratio, h)
    
    G = 6.6743e-11;     %[Nm^2/kg^2] gravitational constant
    
    U = G*M_imp.*M_E./(R_E+R_imp);        %GPE before impact, from r = inf to r = R_E+R_imp
    
    Tf = U*h./((M_E-M_c+M_imp*Si_ratio)*Cp) + T0;   %only incraese in mantle temp

end