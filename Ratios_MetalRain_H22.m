% EL
% July 2022
% Updated 2023-03-06
%
% Iron droplet equilibrium through MO during accretion with given efficiency
%
% Change the geotherm, model type, eq distance, efficiency, composition, file and sheet name
% Comment/uncomment 'writematrix' at the end for output file, as needed

clear;

% PARAMETERS TO CHANGE
geotherm = '3500K';
model = 4;              %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
r_0 = 0.001;            %initial Fe3+/sumFe 
d = 1000;               %[m] eq distance from Rubie+ 2003 (100 to 1000 m)
eff = 0.001;              %[] efficiency factor for droplet equilibrium, 1 = 100%
compSheet = 'Deng20';   %sheet in MoleWeights.xlsx to use for composition
sheetOut = 'H04';           %sheet name to record data
fileOut = 'Rain_H04.xlsx';        % file name


% READ GEOTHERM DATA
data = readmatrix('\db\geotherms.xlsx', 'Sheet', geotherm);
P = data(:,2);
z = data(:,3)*1e3;
T = data(:,4);
r_eq = data(:,7);              %6 = avg, 7 = D20, 8 = A19


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

% pressure at base of MO (assume g = const and base of MO = proportional to
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
J = 2*pi*rho_Si*D*r_d;  %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
% d = 1000;                %[m] eq distance from Rubie+ 2003
dt = d/v_d;             %[s] time interval to fall distance d

% amount of silicate mass equilibriated in time interval dt from iron "rainfall"
dMp_temp = J*dt*N_d*eff; 

r_m_Dt = zeros(1,length(t));        %final r_m over time
r_m_Dt(1) = r_0;
r_m = zeros(1, length(r_eq));       %temp r_m for each iteration
r_m(1) = r_0;
for j = 1:length(dMp_temp)             % earth evolution time
    for i = 2:length(r_eq)        % droplets falling through mantle time
        if P(i) < P_base(j+1)
            dMp = min(dMp_temp(j), Mm*Accr_model(j+1));   %take minimum between equilibrated mass & mantle mass
            % ratio changes as the droplets fall through mantle layers until base of MO
            r_m(i) = (dMp*r_eq(i) + (Mm*Accr_model(j+1)-dMp)*r_m(i-1))/(Mm*Accr_model(j+1));
            idx = i;
        end
    end
    temp = r_m(idx);
    r_m_Dt(j+1) = temp;
    r_m(1) = temp;          % initiate mantle ratio at next iteration with current r
end

% determine logfO2 as dIW for the MO, record surface value.
logfO2vsIW = zeros(1, length(r_m_Dt));
logfO2 = zeros(1, length(r_m_Dt));
for i = 1:length(r_m_Dt)
    [temp_dIW, temp_fO2] = getfO2_H22(P, T, r_m_Dt(i), compSheet);
    logfO2vsIW(i) = temp_dIW(1,1);
    logfO2(i) = temp_fO2(1,1);
end

writematrix(t, fileOut, 'Sheet', sheetOut)
writematrix(r_m_Dt, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)
writematrix(logfO2, fileOut, 'WriteMode', 'append', 'Sheet', sheetOut)


figure(1);
subplot(1,2,1);
hold on
box on
plot(t, r_m_Dt, "LineWidth", 1.5)
yline(r_eq(end), '--', "Label", "r_{eq}="+r_eq(end)+" at base")
yline(mean(r_eq), ':', "Label", "r_{eq} mean ="+mean(r_eq))
ylim([0, r_eq(end)+0.05])
xlabel("Time (Myr)")
ylabel("Fe^{3+}/\SigmaFe Ratio")
hold off

subplot(1,2,2);
hold on
box on
plot(t, logfO2vsIW, "LineWidth", 1.5)
xlabel("Time (Myr)")
ylabel('\DeltaIW (at surface)')

