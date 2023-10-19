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

t_idx = 86;         %test a specific point in time

model = 4;          %accretion model to use, 4 = H04, 5 = N21
%rho_imp = rho_Si*0.68 + rho_Fe*0.32;            % weighted average for density of impactor
rho_imp = 5000;     %[kg/m^3] approximation based on weighted average (above)
M_E_0 = 5.97e24;    %[kg] present day mass of Earth
M_c_0 = 1.88e24;    %[kg] present day mass of core
R = 8.314;          %[J/mol*K] gas const 

%note: only have this set up to do the Pmo method

d_mo_factor = 1;    %choose depth of melting, where 0 = surf and 1 = CMB

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

%read spreadsheet of dV as function of P for given constant T
data = readmatrix('/db/dVvsP.xlsx');
dV_test = data(:,3:10)*1e-6;        %[m^3/mol] is in cm^3 in spreadsheet
P_test = data(:,2)*1e9;
T_test = [2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500];

%reshape T to match dV data P and T
T(isnan(T))=[];
idx_length = length(T);
T_f = interp1(P(1:idx_length), T, P_test);

PV = zeros(length(P_test),1);

for i = 2:length(P_test)
    if isfinite(T_f(i))
        dV_q = zeros(i,1);
        for k = 1:i                 %query for dV at constant T(i), index k matches row in P
            dV_q(k) = interp1(T_test, dV_test(k,:),T_f(i), 'pchip');
            %can't use interp2 because can't extrapolate(?)
            %you will get interp1 errors on ignoring NaN, that's okay with 'pchip' method
        end
        int = trapz(P_test(1:i),dV_q);
        PV(i) = int/(R*T_f(i));           %unitless quantity
    else
        PV(i) = NaN;
    end
end

T_f2 = T_f;
T_f2(isnan(T_f2))=[];
idx_length = length(T_f2);
PV_f = interp1(T_f2, PV(1:idx_length), T, 'linear', 'extrap');

% test with function
%PV_f2 = calcPV(T,P);

disp(Tp(t_idx))

% data in PVcalc is based on geotherms from M&K2019b Fig 1
data = readmatrix('PVcalc.xlsx');
P_check = data(:,1);
PV_check = data(:,2:end);


figure(1);
hold on
box on
plot(PV_f,P(1:length(PV_f))/1e9, "LineWidth",2)
%plot(PV_f2,P(1:length(PV_f2))/1e9, 'r--', "LineWidth", 2)
plot(PV_check,P_check)