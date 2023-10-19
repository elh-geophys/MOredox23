% EL
% July 2022
% Updated 2023-08-31
%
% Iron droplet equilibrium through MO during accretion with given efficiency
%
% Change the geotherm, model type, eq distance, efficiency, composition, file and sheet name
% Comment/uncomment 'writematrix' at the end for output file, as needed


% TO FIX:  1) update geotherm calc with new accretion modeling
%          2) query PV calc given geotherm
%          3) calc Fe ratio (with new PV calc and varying IW based on accretion)
%          4) test and test metal rain
%          5) update the MC metal rain accordingly

clear;

% PARAMETERS TO CHANGE
geotherm = '3500K';
model = 5;                      %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
r_0 = 0.001;                    %initial Fe3+/sumFe 
eff = 0.5;                      %[] efficiency factor for droplet equilibrium, 1 = 100%
r_eqSheet = 'N21_3500K';        %sheet for r_eq(t) data
compSheet = 'Deng20';           %sheet in MoleWeights.xlsx to use for composition
sheetOut = 'H04';                %sheet name to record data
fileOut = 'Rain.xlsx';        % file name
write = 0;                    %1 or 0, to write to file

% READ GEOTHERM DATA
data = readmatrix('\db\geotherms.xlsx', 'Sheet', geotherm);
P = data(:,2);
z = data(:,3)*1e3;
T = data(:,4);

% READ R_eq DATA
data = readmatrix('\db\r_eq_D20.xlsx', 'Sheet', r_eqSheet);
r_eq = data(2:end,2:end);              


% DETERMINATION OF FE RATIO & IW DURING ACCRETION
Mm = 4e24;        %[kg] mass of whole mantle
R_E_0 = 6371e3;   %[m] radius of Earth
rho_E = 5500;     %[kg/m^3] mean density of Earth
rho_Si = 3750;    %[kg/m^3] silicate melt density

% First 100 Myr accretion history of Earth
[t,Accr_model] = getAccrModel(model);

% estimate time evolution of iron on Earth
M_earth = 5.97e24;              % [kg] mass of Earth
Fe_earth = 0.321;               % [] weight % of iron
M_Fe_earth = M_earth*Fe_earth;  % [kg] estimated mass of iron on Earth
M_Fe_accr = M_Fe_earth*Accr_model;   % [kg] time evolve of mass of iron

% radius of Earth over time
R_E = (3*Accr_model*M_earth/(4*pi*rho_E)).^(1/3);

% pressure at base of MO (assume base of MO = proportional to
% present day mantle depth)
z_base = z(end)*R_E/R_E_0;
P_base = interp1(z,P,z_base);

% iron droplets
r_d = 0.5e-2;                   %[m] radius of droplets (~1cm diameter)
v_d = 0.5;                      %[m/s] percolation velocity
rho_Fe = 7800;                  %[kg/m^3] iron melt density
m_d = rho_Fe * 4/3*pi*r_d^3;    %[kg] mass of a single droplet
dM_Fe = diff(M_Fe_accr,1,2);    %[kg] iron mass accreted per time interval
N_d = dM_Fe/m_d;                % # of droplets per time interval

D = 1e-7;               %[m^2/s] chemical diffusity

% Using eq distance from Rubie+ 2003
d = 1000;               %[m] eq distance from Rubie+ 2003
dt = d/v_d;             %[s] time interval to fall distance d

% Ulvrova+2011, diffusion boundary for "drop" with high Re and low viscosity ratio
a = 0.79;
Pe = r_d*v_d/D;     %Peclet #
k = 0.1;            %partition coefficient
Sh = a*Pe^(0.5);    %Shmidt #
R_D = 1;            %diffusivity ratio ~ 1  (may be 1-10?)
h = r_d/(a*Pe^(0.5));                   %diffusion boundary
tau = Pe/3 * ( k/Sh + 1/(10*R_D)) * r_d/v_d;      %time scale for equilibrium

%mass flow rate J = A*v*rho
J_old = 2*pi*rho_Si*D*r_d;          %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
J = 2*pi*r_d*h*v_d*rho_Si/2;    %[kg/s] "droplet ring" mass flow rate, divide by 2 for diffusion gradient

% amount of silicate mass equilibriated in time interval dt (or tau) from iron "rainfall"
dMp_temp_old = J_old*dt*N_d*eff;                                           %cylindrical tail or ring model
dMp_temp = J*tau*N_d*eff;                                           %spherical shell, with time

r_m_Dt = zeros(1,length(t));        %final r_m over time
r_m_Dt(1) = r_0;
r_m = zeros(1, length(P));       %temp r_m for each iteration
r_m(1) = r_0;
for j = 1:500             % earth evolution time
    for i = 2:length(P)        % droplets falling through mantle time
        if P(i) < P_base(j+1)
            dMp = min(dMp_temp(j), Mm*Accr_model(j+1));   %take minimum between equilibrated mass & mantle mass
            % ratio changes as the droplets fall through mantle layers until base of MO
            r_m(i) = (dMp*r_eq(i,j) + (Mm*Accr_model(j+1)-dMp)*r_m(i-1))/(Mm*Accr_model(j+1));
            idx = i;
        end
    end
        figure(4);
        hold on
        plot(r_eq(:,j), 'k--')
        plot(r_m, 'm')
        
    temp = r_m(idx);
    r_m_Dt(j+1) = temp;
    r_m(1) = temp;          % initiate mantle ratio at next iteration with current r
end

% determine logfO2 as dIW for the MO, record surface value.
% logfO2vsIW = zeros(1, length(r_m_Dt));
% logfO2 = zeros(1, length(r_m_Dt));
% for i = 1:length(r_m_Dt)
%     [temp_dIW, temp_fO2] = getfO2_H22(P, T, r_m_Dt(i), compSheet);
%     logfO2vsIW(i) = temp_dIW(1,1);
%     logfO2(i) = temp_fO2(1,1);
% end

if write == 1
    writematrix(t, fileOut, 'Sheet', sheetOut)
    writematrix(r_m_Dt, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)
    %writematrix(logfO2, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)
end

figure(1);
subplot(1,2,1);
hold on
box on
plot(t, r_m_Dt, "LineWidth", 1.5)
yline(r_eq(end), '--', "Label", "r_{eq}="+r_eq(end)+" at base")
yline(mean(r_eq(:,end)), ':', "Label", "r_{eq} mean ="+mean(r_eq(:,end)))
ylim([0, r_eq(end)+0.05])
xlabel("Time (Myr)")
ylabel("Fe^{3+}/\SigmaFe Ratio")
hold off

subplot(1,2,2);
plot(t(2:end), dMp_temp)
xlabel("Time(Myr)")
ylabel("Mass Equilibrated (kg)")
% hold on
% box on
% plot(t, logfO2vsIW, "LineWidth", 1.5)
% xlabel("Time (Myr)")
% ylabel('\DeltaIW (at surface)')

figure(2);
hold on
plot(dMp_temp)
plot(Mm*Accr_model)
xlabel("Time Myr")
ylabel("Mass (kg)")
legend("Eq'd Si", "Mantle Mass")

