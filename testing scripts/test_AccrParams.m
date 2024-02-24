% EL
% August 2023

% Modeling accretion of growing Earth
% Pressure-depth dependence and geotherm

clear;

model = 4;          %accretion model to use, 4 = H04, 5 = N21
rho_Si = 3750;      %[kg/m^3] silicate melt density
rho_Fe = 7800;      %[kg/m^3] iron density
%rho_imp = rho_Si*0.68 + rho_Fe*0.32;            % weighted average for density of impactor
rho_imp = 5000;     %[kg/m^3] approximation based on weighted average (above)
M_E_0 = 5.97e24;    %[kg] present day mass of Earth
M_c_0 = 1.88e24;    %[kg] present day mass of core
G = 6.6743e-11;     %[Nm^2/kg^2] gravitational constant
Cp = 1e3;           %[J/kgK] specific heat 
T0 = 1613;          %[K] assume T starts at rheological transition 1613K Tp

h = 0.2;            %energy contribution factor for U
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
R_imp = (3/(4*pi)*M_imp/rho_imp).^(1/3);    %radius of impactor 

% So plan... use average core and mantle radius to determine the possible
% MO depth, BUT use all of the mantle mass for chemical mixing after each
% iterative step in oxidation code.

% JK says to use average core mass from before and after each impact 

%T&S Geodyanm Eqn 2.73, Pressure as a function radius from center with
%varying size of planet
%P_cmb = 4/3*pi*G*rho_m.*R_c.^3.*(rho_c-rho_m).*(1./R_c-1./R_E) + 2/3*pi*G*rho_m.^2.*(R_E.^2-R_c.^2); 

%test function
P_cmb = zeros(1,length(R_c));
for i = 1:length(R_c)
    P_cmb(i) = get2LayerP(rho_m(i), rho_c(i), R_E(i), R_c(i), R_c(i));
end
P_cmb(isnan(P_cmb))=0;

% choose impactor MO depth
P_mo = P_cmb*d_mo_factor;

% for geotherm, two variations: 1) Tp based on MO depth, 2) change in temp based on impact PE

% variation 1: Tp based on MO depth
% 
% data = readmatrix('\db\geotherms_combo.xlsx');
% P_data = data(:,1);
% Tliq = data(:,4);
% T_data = data(:,6:10);
% 
% Tliq_cross = zeros(1,5);
% for i = 1:5
%     [~,idx] = min(abs(T_data(:,i)-Tliq));
%     Tliq_cross(i) = T_data(idx,i);
% end
% 
% T_test = [2000, 2500, 3000, 3500, 4000];
% P_test = P_mo/1e9;                          %testing MO pressure
% 
% Tliq_test = interp1(P_data,Tliq,P_test, 'pchip', 'extrap');
% Tliq_test_cmb = interp1(P_data,Tliq,P_cmb/1e9, 'pchip', 'extrap');
% 
% Tp = interp1(Tliq_cross,T_test,Tliq_test, 'pchip', 'extrap');
% Tp_cmb = interp1(Tliq_cross,T_test,Tliq_test_cmb, 'pchip', 'extrap');

Tp = getMOTp(P_mo/1e9);
Tp_cmb = getMOTp(P_cmb/1e9);

P_max = round(P_mo,-9);
Tad = zeros((P_max(end)/1e9*2)+1,length(P_max));      %set up Tad

%determine adiabat, from code by J.Korenaga
%Tad = zeros(length(P),1);
P = ((0:1e9:(P_max(1)*2))/2)';           %to go by 0.5 GPa
if size(P)==1
    P = [0, 0.5e9];
end

% Tad(1,1) = Tp(1);
% dp = P(2)-P(1);
% for j=2:length(P)                        %for the first one
%   alpha = 3.622e-5*exp(-2.377e-5*Tad(j-1,1)-0.0106*P(j-1)/1e9);
%   cp =  627+0.441*Tad(j-1,1)-0.211*P(j-1)/1e9;
%   rho = 2870-0.082*Tad(j-1,1)+162*(P(j-1)/1e9)^0.58;
%   dTdp = alpha*Tad(j-1,1)/(rho*cp);      
%   dT = dTdp*dp;
%   Tad0 = Tad(j-1,1)+dT;
%   Tad(j,1) = Tad0;
% end

Tad_temp = getMOAdiabat(Tp(1),P);
Tad_empty = zeros(size(Tad,1)-size(Tad_temp,1),1);
Tad_temp2 = cat(1, Tad_temp, Tad_empty);
Tad = Tad_temp2;


for i = 2:length(P_max)
    if M_imp(i-1) > 0
        P = ((0:1e9:(P_max(i)*2))/2)';           %to go by 0.5 GPa
        
        Tad_temp = getMOAdiabat(Tp(i),P);
        Tad_empty = zeros(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);
        
    else
        
        Tad_temp = Tad(:,i-1);
        Tad_empty = zeros(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);
    end
   
end

% variation 2: Tp based on impact PE

M_imp_temp = cat(2,M_imp,[0]);      %so it matches length of M_E and R_E
R_imp_temp = cat(2,R_imp,[0]);
Tf_test_temp = calcTforU2Q(T0,Cp,M_E,M_c,M_imp_temp,R_E,R_imp_temp, 0.679, h);

Tf_test = zeros(1,length(M_E));
Tf_test(1) = T0;
for i = 2:length(Tf_test_temp)
    if Tf_test_temp(i-1) > Tp(i)
        Tf_test(i) = Tp(i);
    else
        Tf_test(i) = Tf_test_temp(i-1);
    end
end

Tf = Tf_test;

P = ((0:1e9:(P_max(1)*2))/2)';           %to go by 0.5 GPa
if size(P)==1
    P = [0, 0.5e9];
end

TUad = zeros((P_max(end)/1e9*2)+1,length(P_max));      %set up Tad

TUad_temp = getMOAdiabat(Tf(1),P);
TUad_empty = zeros(size(TUad,1)-size(TUad_temp,1),1);
TUad_temp2 = cat(1, TUad_temp, TUad_empty);
TUad = TUad_temp2;


for i = 2:length(P_max)
    if M_imp(i-1) > 0
        P = ((0:1e9:(P_max(i)*2))/2)';           %to go by 0.5 GPa
    
        TUad_temp = getMOAdiabat(Tf(i),P);
        TUad_empty = zeros(size(TUad,1)-size(TUad_temp,1),1);
        TUad_temp2 = cat(1, TUad_temp, TUad_empty);
        TUad = cat(2,TUad,TUad_temp2);
        
    else        
        TUad_temp = TUad(:,i-1);
        TUad_empty = zeros(size(TUad,1)-size(TUad_temp,1),1);
        TUad_temp2 = cat(1, TUad_temp, TUad_empty);
        TUad = cat(2,TUad,TUad_temp2);
    end
   
end



figure(1);
hold on
box on
for i = 2:length(M_E)-1
    if M_imp(i) > 0
        plot(t(i+1),Tp(i+1), 'bo', 'MarkerSize', 8)
        plot(t(i+1),Tf(i+1), 'rx', 'MarkerSize', 8)
    end
end
xlabel('Time (Myr)')
ylabel('Potential Temperature (K)')
ylim([1500 4500])
xlim([0 60])
legend('Based on CMB @ Liquidus', 'Based on U -> Heat', 'Location', 'southeast')
title(['Model ', num2str(model), ' with f=', num2str(d_mo_factor)])

figure(2);
hold on
box on
for i = 2:length(M_E)-1
    if M_imp(i) > 0
        plot(Tad(:,i+1),P/1e9, 'b-')
        plot(TUad(:,i+1),P/1e9, 'r-')
    end
end
xlabel('Temp(K)')
ylabel('Pressure (GPa)')
ylim([0 135])
xlim([1000 6000])
legend('Based on CMB @ Liquidus', 'Based on U -> Heat', 'Location', 'southeast')
title(['Model ', num2str(model), ' with f=', num2str(d_mo_factor)])

