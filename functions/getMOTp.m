% EL
% Aug 31, 2023
%
% Function to get Tp in magma ocean for where the base of the magma ocean
% meets liquidus.
%
% REQUIRED: geotherms_combo.xlsx, from M&K2019b/K2023
%
% INPUTS:    P_mo_base     [Pa] the pressure at the MO base, scalar or 1D array
%            data          data from geotherms_combo.xlsx, created by code from J.Korenaga for magma ocean adiabats at different Tp
%
% OUTPUTS:   Tp            [K] potential temperature, scalar or 1D array
%

function [Tp] = getMOTp(P_mo_base, data)
    %data = readmatrix('\db\geotherms_combo.xlsx');
    P_data = data(:,1);
    Tliq = data(:,4);
    T_data = data(:,6:10);

    Tliq_cross = zeros(1,5);
    for i = 1:5
        [~,idx] = min(abs(T_data(:,i)-Tliq));
        Tliq_cross(i) = T_data(idx,i);
    end

    T_test = [2000, 2500, 3000, 3500, 4000];
    P_test = P_mo_base;                          %testing MO pressure

    Tliq_test = interp1(P_data,Tliq,P_test, 'pchip', 'extrap');    %the liquidus pt that matches P_mo

    Tp = interp1(Tliq_cross,T_test,Tliq_test, 'pchip', 'extrap');  %interpolated Tp
end

