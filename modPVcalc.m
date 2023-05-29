% EL
% modified PV calc

% So. Here's the catch. The integration term is dV @ that given
% temperature for the *WHOLE* P0 to P.  It is *NOT* dV along geotherm from
% P0 to P... as we assumed.  So how do we fix this?  Well, Deng did the dV
% calc along the geotherm for *each* individual temperature, then stored PV
% into excel and used that for all further fO2 calcs.  There's a sub
% routine called PV_cal_geo which does this.  However, my geotherm will not
% for whatever weird reason actually come up with reasonable results past
% 4000K using this routine.  It's unclear why.  As a compromise, I've
% established dV along constant temps using Deng Python code (as shown in 
% Fig 1 of their paper). This is stored in dVvsP.xlsx.

% The plan here is to pull dV vs T from excel sheet for each point P,
% interpolate/extrapolate to the given T, and output dV to integrate.

clear;

data = readmatrix('/db/geotherms_combo.xlsx');
P = data(:,1)*1e9;
T = data(:,11);         %choose geotherm

data = readmatrix('/db/dVvsP.xlsx');
dV_test = data(:,3:10)*1e-6;
T_test = [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500];

figure(4);
plot(P/1e9, dV_test);

R = 8.314;      %gas const

PV = zeros(length(P),1);

for i = 2:length(P)
    dV_q = zeros(i,1);
    for j = 1:i
        dV_q(j) = interp1(T_test, dV_test(j,:),T(i),'pchip');
    end
    
%     disp(dV_q);
%     disp(dV(1:i));
%     disp(T(i));
%     pause(0.5);
     
%     figure(1);
%     plot(T_test,dV_test(i,:),T(i),dV_q(i),'ro');
%     title("End point, P = "+P(i)/1e9)
%     
%     r = round(i/2);
%     figure(2);
%     plot(T_test,dV_test(r,:),T(i),dV_q(r),'ro');
%     title("Mid point, P = "+P(r)/1e9)
%     
    int = trapz(P(1:i),dV_q);
    PV(i) = int/(R*T(i));

end


% writematrix(P/1e9, 'PVcalc.xlsx')
% writematrix(PV, 'PVcalc.xlsx', 'Range', 'G1') 


figure(3);
hold on
box on
plot(PV,P,'r-')
hold off