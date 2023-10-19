% EL
% Sept 19, 2023

% Determine PV integration as a function of P and T to generate an output
% file.
% Testing some different PV interpolation routines also.

% to turn off and on warnings from calcPV.m:
%  warnID = 'MATLAB:interp1:NaNstrip'; warning('off', warnID)
%  warning('on')

clear;

fileOut = 'PVcalc.xlsx';
write = 1;

Tmin = 2000;
Tmax = 5500;

T = Tmin:100:Tmax;

R = 8.314;          %[J/mol*K] gas const

dVvsP_data = readmatrix('/db/dVvsP.xlsx');
Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');

%read spreadsheet of dV as function of P for given constant T
dV_test = dVvsP_data(:,3:10)*1e-6;        %[m^3/mol] is in cm^3 in spreadsheet
P_test = dVvsP_data(:,2)*1e9;
T_test = [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500];

PV_data = zeros(length(P_test),length(T));

for j = 1:length(T)
    for i = 2:length(P_test)            %P_test and T_f match arrays (temp as function of pressure)
        dV_q = zeros(i,1);
        for k = 1:i
            if P_test(k) < 0.5e9 && T(j) <= 4000
                dV_q(k) = interp1(T_test(1:5), dV_test(k,1:5),T(j),'pchip');
            elseif P_test(k) < 10e9 && T(j) <= 4500
                dV_q(k) = interp1(T_test(1:6), dV_test(k,1:6),T(j),'pchip');
            else
                % withNaN = sum(isnan(dV_test(k,:)));
                % if withNaN > 0
                % disp(['interp with NaN at ', num2str(P_test(k)/1e9)])
                % end
                dV_q(k) = interp1(T_test, dV_test(k,:),T(j),'pchip');         %this is the original
                %disp(dV_test(k,:))
            end
        end

        int = trapz(P_test(1:i),dV_q);
        PV_data(i,j) = int/(R*T(j));           %unitless quantity

    end
end

if write == 1
    writematrix(P_test/1e9, fileOut, 'Range', 'A2');
    writematrix(T, fileOut, 'Range', 'B1');
    writematrix(PV_data, fileOut, 'Range', 'B2');
end

Tp = 3500;

dP = 0.5e9;
P = (0:dP:126e9)';


Tad = getMOAdiabat(Tp,P,Adiabat_data);

% DETERMINE PV INTEGRATION AS INT(PV)/RT
PV_byfunc = calcPV(Tad,P, dVvsP_data);

PV_byinterp2D = interp2(T, P_test, PV_data, Tad, P);

Tad_test = interp1(P, Tad, P_test);
Tad_test(isnan(Tad_test))=[];
PV_byinterp1D_temp = zeros(length(Tad_test),1);
for i = 1:length(Tad_test)
    PV_byinterp1D_temp(i) = interp1(T, PV_data(i,:), Tad_test(i), 'pchip');
end

PV_byinterp1D = interp1(Tad_test, PV_byinterp1D_temp, Tad, 'linear', 'extrap');

PercDiff1 = (PV_byfunc - PV_byinterp2D)./PV_byfunc * 100;
PercDiff2 = (PV_byfunc - PV_byinterp1D)./PV_byfunc * 100;

figure(1);
hold on
plot(P, PercDiff1)
plot(P, PercDiff2)

figure(2);
hold on
plot(P, PV_byfunc)
plot(P, PV_byinterp2D)
plot(P, PV_byinterp1D)

