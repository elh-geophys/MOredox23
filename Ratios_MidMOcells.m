% EL
% June/July 2022
% Updated 2023-03-13
%
% Determine the evolution of Fe ratio and logfO2 by large convection
% cell mixing (mid-mantle MO version)

clear;

Tp = '3500K';

data = readmatrix('\db\geotherms.xlsx', 'Sheet', 'PREM');
z = data(:,3)*1e3;

data = readmatrix('\db\geotherms.xlsx', 'Sheet', Tp);
P = data(:,2);
T = data(:,4);
r_eq = data(:,7);           %6 = avg, 7 = D20, 8 = A19

r_0 = 0.001;

% PARAMETERS
D = 1e-7;       %[m^2/s] chemical diffusivity
n = 20;         %[] number of convection cells
Tm = 3500;      %[K] mantle potential temperature (check 'sheet', line 15)
Pbase = 100;    %[GPa] base of MO (25, 50, 80, 100, 135)
r_eq_base = interp1(P,r_eq,Pbase);      %ratio at base of MO

% set the time interval
dt = 1e4;       %[yr]
tmax = 5e6;
N = tmax/dt+1;
t_yr = linspace(0, tmax, N);
t = t_yr*365*24*60*60;          %[s]

% CONSTANTS
Rc = 3400e3;            %[m] core radius
Mm = 4e24;              %[kg] mass of whole mantle
alpha = 5e-5;           %[1/K] thermal expansivity
g = 9.8;                %[m/s^2] gravitational acceleration
rho = 3750;             %[kg/m^3] silicate melt density
Cp = 1e3;               %[J/kgK] heat capacity
visc = 0.1;             %[Pa s] viscosity
k = 2;                  %[W/mK] thermal conductivity
kappa = k/(rho*Cp);     %[m^2/s] thermal diffusivity
k_visc = visc/rho;      %[m^2/s] kinematic viscosity
C = 2*pi*Rc;            %[m] circumference of core (1D)   

% surface temperature based on mantle potential temp
c0 = 546;           %from parameterization M&K2022
c1 = 0.63;
Ts = c0 + c1*Tm;    %[K] surface temp

L = interp1(P,z,Pbase); %[m] depth of MO

Ra = alpha*g*(Tm-Ts)*L^3/(kappa*k_visc);    %Ra number

F = 0.089*k*(Tm-Ts)*Ra^(1/3)/L;             %surface flux, Solomatov 2015

v = 0.6*(alpha*g*L*F/(rho*Cp))^(1/3);       %soft-conv velocity, Solomatov 2015

R = Rc + (2890e3-L);           %[m] distance from center of Earth of base of MO
l = sqrt(4*pi*R^2/n);          %[m] length scale of area broken up into n pieces
J = 4*pi*rho*R^2*sqrt(D*v/l);  %[kg/s] mass flux, simplified from J = rho*Vol/tau

Vs = 4/3*pi*(R^3 - Rc^3);             %[m^3] volume of solid mantle shell
Vmo = 4/3*pi*((Rc+z(end))^3 - R^3);   %[m^3] volume of MO
dMp = J*(t(2)-t(1));                  %[kg] mass processed in time interval dt
Mmo = Mm * Vmo/(Vmo+Vs);              %[kg] mass of MO

r = zeros(1,length(t));
r(1) = r_0;
for i = 2:length(t)
    r(i) = (dMp*r_eq_base + (Mmo-dMp)*r(i-1))/Mmo;
end

% %       0yr,  10kyr, 100kyr,500kyr,1Myr,  5Myr
% r_div = [r(1), r(2), r(11), r(51), r(101), r(501)];
% % IW as function of r and T
% logfO2vsIW = zeros(length(T), length(r_div));
% for i = 1:length(r_div)
%     logfO2vsIW(:,i) = fO2_EL(P, T, dV, r_div(i), G_flag, PV_method);
% end


% fO2 as function of time, record only surface value
logfO2_surf = zeros(1, length(r));
for i = 1:length(r)
    [~, logfO2] = getfO2_H22(P,T,r(i),'Deng20');
    logfO2_surf(i) = logfO2(1);
end

writematrix(t_yr, 'MidMOcells.xlsx', 'Sheet', Tp)
writematrix(r, 'MidMOcells.xlsx', 'WriteMode', 'append', 'Sheet', Tp)
writematrix(logfO2_surf, 'MidMOcells.xlsx', 'WriteMode', 'append', 'Sheet', Tp)

figure(1);
hold on
plot(t_yr,r)
yline(r_eq_base,'--')
ylim([0 0.3])
xlabel("t (yr)")
ylabel("ratio")
hold off

figure(2);
hold on
plot(t_yr, logfO2_surf, "LineWidth", 1.5)
xlabel("t (yr)")
ylabel('Redox (logf_{O_2})')

% figure(3);
% hold on
% box on
% newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
% colororder(newcolors)
% plot(logfO2vsIW,P, 'LineWidth', 1.5)
% ylabel('Pressure (GPa)')
% xlabel('Redox (logf_{O_2}– IW)')
% xlim([-20 4])
% ylim([0 Pbase])
% xticks(-20:2:6)
% legend(["0 yr", "10 kyr", "100 kyr", "500 kyr",  "1 Myr", "5 Myr"], 'Location', 'northeastoutside')
% ax = gca;
% ax.YDir = 'reverse';