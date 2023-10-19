% EL
% Sept 5, 2023
%
% Function to calculate PV integral (work term). The integration of PdV goes from P0 to
% P, as a constant temperature. Using Python code given to me by J. Deng, I
% have generated an excel file with dV as a function of P at constant
% temperatures (dVvsP.xlxs). Using this, we can interpolate/extrapolate
% values for any given P and T to calculate int(PdV). This function takes
% in a geotherm (T as function of P) and calculates PV. See test script
% modPVcalc.m for more information.
%
% INPUTS:   T       [K] geotherm temperature, 1D array
%           P       [Pa] geotherm pressure, 1D array
%
% OUTPUTS:  PV_f    [] integration term as int(PdV)/RT
%
% This process was streamlined in /testing scripts/PVintegration.m to
% generate an output file.
%

function [PV_f] = calcPV_old(T,P,data)

    R = 8.314;          %[J/mol*K] gas const

    %read spreadsheet of dV as function of P for given constant T
    %data = readmatrix('/db/dVvsP.xlsx');
    dV_test = data(:,3:10)*1e-6;        %[m^3/mol] is in cm^3 in spreadsheet
    P_test = data(:,2)*1e9;
    T_test = [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500];

    % shape T to match test data
    T(isnan(T))=[];
    idx_length = length(T);
    T_f = interp1(P(1:idx_length), T, P_test);

    PV = zeros(length(P_test),1);

    for i = 2:length(P_test)            %P_test and T_f match arrays (temp as function of pressure)
        if isfinite(T_f(i))
            dV_q = zeros(i,1);
            for k = 1:i
                if P_test(k) < 0.5e9 && T_f(i) <= 4000
                    %can't use interp2 because can't extrapolate(?)
                    dV_q(k) = interp1(T_test(1:5), dV_test(k,1:5),T_f(i),'pchip');
                elseif P_test(k) < 10e9 && T_f(i) <= 4500
                    dV_q(k) = interp1(T_test(1:6), dV_test(k,1:6),T_f(i),'pchip');
                else
                    % withNaN = sum(isnan(dV_test(k,:)));
                    % if withNaN > 0
                    % disp(['interp with NaN at ', num2str(P_test(k)/1e9)])
                    % end
                    dV_q(k) = interp1(T_test, dV_test(k,:),T_f(i),'pchip');         %this is the original
                    %disp(dV_test(k,:))
                end
            end

            int = trapz(P_test(1:i),dV_q);
            PV(i) = int/(R*T_f(i));           %unitless quantity
                      
        else
            PV(i) = NaN;        %save on computation time
        end
    end

    % reshape to match inputs
    T_f(isnan(T_f))=[];
    idx_length = length(T_f);
    PV_temp = interp1(T_f, PV(1:idx_length), T, 'linear', 'extrap');
    
    PV_nan = NaN(length(P)-length(PV_temp),1);
    PV_f = cat(1,PV_temp,PV_nan);                   %so retains NaNs at output to match length of T array input (if NaNs to start)

end
