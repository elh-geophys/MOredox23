% EL
% Sept 5, 2023
%
% Function to calculate PV integral (work term). The integration of PdV goes from P0 to
% P, as a constant temperature. Using Python code given to me by J. Deng, I
% have generated an excel file with dV as a function of P at constant
% temperatures (dVvsP.xlxs). Using this, we can interpolate/extrapolate
% values for any given P and T to calculate int(PdV). This is done by PVintegration.m
% % This function takes in a geotherm (T as function of P) and determine
% int(PV)/RT from interp/extrap of data from PV_calc.xlsx.
%
% INPUTS:   T       [K] geotherm temperature, 1D array
%           P       [Pa] geotherm pressure, 1D array
%           PVdata  [] taken from PV_calc.xlsx
%
% OUTPUTS:  PV      [] integration term as int(PdV)/RT
%

function [PV] = calcPV(T,P, PVdata)

    P_test = PVdata(2:end,1)*1e9;
    T_test = PVdata(1,2:end);
    PV_test = PVdata(2:end,2:end);        %note: this is int(PdV)/RT, unitless quantity

    Tad_test = interp1(P, T, P_test);
    Tad_test(isnan(Tad_test))=[];
    
    PV_temp = zeros(length(Tad_test),1);
    for i = 1:length(Tad_test)
        PV_temp(i) = interp1(T_test, PV_test(i,:), Tad_test(i), 'pchip');
    end

    PV = interp1(Tad_test, PV_temp, T, 'linear', 'extrap');
    
end
