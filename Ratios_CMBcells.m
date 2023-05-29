% EL
% June/July 2022
% Updated 2023-03-13
%
% Determine the evolution of Fe ratio and logfO2 by large convection
% cell mixing.

clear;

Tp = '4500K';

data = readmatrix('\db\geotherms.xlsx', 'Sheet', Tp);
P = data(:,2);
T = data(:,4);
r_eq = data(end,7);           %6 = avg, 7 = D20, 8 = A19

r_0 = 0.001;

% PARAMETERS
D = 1e-7;       %[m^2/s] chemical diffusivity
n = 20;         %[] number of convection cells
Tm = 4500;      %[K] mantle potential temperature (check 'sheet', line 9)


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
L = 2880e3;             %[m] depth of mantle
rho = 3750;             %[kg/m^3] silicate melt density
Cp = 1e3;               %[J/kgK] heat capacity
visc = 0.1;             %[Pa s] silicate melt viscosity
k = 2;                  %[W/mK] thermal conductivity
kappa = k/(rho*Cp);     %[m^2/s] thermal diffusivity
k_visc = visc/rho;      %[m^2/s] kinematic viscosity
C = 2*pi*Rc;            %[m] circumference of core (1D)   

% surface temperature based on mantle potential temp
c0 = 546;           %from parameterization M&K2022
c1 = 0.63;
Ts = c0 + c1*Tm;    %[K] surface temp

Ra = alpha*g*(Tm-Ts)*L^3/(kappa*k_visc);    %Ra number

F = 0.089*k*(Tm-Ts)*Ra^(1/3)/L;             %surface flux, Solomatov 2015

v = 0.6*(alpha*g*L*F/(rho*Cp))^(1/3);       %soft-conv velocity, Solomatov 2015

l = sqrt(4*pi*Rc^2/n);          %[m] length scale of area broken up into n pieces
J = 4*pi*rho*Rc^2*sqrt(D*v/l);  %[kg/s] mass flux, simplified from J = rho*Vol/tau

dMp = J*(t(2)-t(1));            %[kg] mass processed in time interval dt

r = zeros(1,length(t));
r(1) = r_0;
for i = 2:length(t)
    r(i) = (dMp*r_eq + (Mm-dMp)*r(i-1))/Mm;
end

% fO2 as function of time, record only surface value
logfO2_surf = zeros(1, length(r));
for i = 1:length(r)
    [~, logfO2] = getfO2_H22(P,T,r(i),'Deng20');
    logfO2_surf(i) = logfO2(1);
end

writematrix(t_yr, 'CMBcells.xlsx', 'Sheet', Tp)
writematrix(r, 'CMBcells.xlsx', 'WriteMode', 'append', 'Sheet', Tp)
writematrix(logfO2_surf, 'CMBcells.xlsx', 'WriteMode', 'append', 'Sheet', Tp)

figure(1);
hold on
plot(t_yr,r)
yline(r_eq,'--')
ylim([0 0.3])
xlabel("t (yr)")
ylabel("ratio")
hold off

figure(2);
hold on
plot(t_yr, logfO2_surf, "LineWidth", 1.5)
xlabel("t (yr)")
ylabel('log(f_{O2})')
