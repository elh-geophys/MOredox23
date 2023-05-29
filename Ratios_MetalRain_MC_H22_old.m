% EL
% August 2022
% Updated 2023-02-20
%
% Produce Monte Carlo results with varying efficiency and depth for each
% impact during accretion.
% Uses getRainRatio.m

clear;

data = readmatrix('\db\geotherms.xlsx', 'Sheet', '3500K');
P = data(:,2);
z = data(:,3)*1e3;
T = data(:,4);
r_eq = data(:,7);           %6 = avg, 7 = D20, 8 = A19

r_0 = 0.001;
compSheet = 'Deng20';   %composition to use
model = 4;          %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21

sheet = 'H04_r';                %sheet name to record data
sheetIW = 'H04_dIW';
sheetfO2 = 'H04_fO2';
file = 'Rain_MC.xlsx';        % file name

% CONSTANTS
Mm = 4e24;        %[kg] mass of whole mantle
R_E_0 = 6371e3;   %[m] radius of Earth
rho_E = 5500;     %[kg/m^3] mean density of Earth
rho_Si = 3750;    %[kg/m^3] silicate density
D = 1e-7;         %[m^2/s] chemical diffusity
d = 1000;         %[m] eq distance from Rubie+ 2003 (100 to 1000 m)

% First 100 Myr accretion history of Earth
[t, Accr_model] = getAccrModel(model);

% estimate time evolution of iron on Earth
M_earth = 5.97e24;                  % [kg] mass of Earth
Fe_earth = 0.321;                   % [] weight % of iron
M_Fe_earth = M_earth*Fe_earth;      % [kg] estimated mass of iron on Earth
M_Fe_accr = M_Fe_earth*Accr_model;  % [kg] time evolve of mass of iron

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

J = 2*pi*rho_Si*D*r_d;  %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
dt = d/v_d;             %[s] time interval to fall equilibrium distance d
dMp_temp = J*dt*N_d;

% INPUTS FOR FUNCTION
% P, z, z_base, R_E, rho_Si, Mm, Accr_model, r_eq, dMp_temp, Mm, r_0, t

r_m_Dt = zeros(1000,length(t));
for i = 1:length(t)
	r_m_Dt(i,:) = getRainRatio(P,z,r_0,r_eq,t,Accr_model,z_base,R_E,rho_Si,Mm,dMp_temp);
end

% Determine the average and mid-50% and mid-90% 
r_avg = median(r_m_Dt,1);
r_prc = prctile(r_m_Dt,[5 25 75 95],1);
r_5p =  r_prc(1,:);
r_25p = r_prc(2,:);
r_75p = r_prc(3,:);
r_95p = r_prc(4,:);

r = [r_avg; r_5p; r_25p; r_75p; r_95p];

% determine logfO2 as dIW for the MO, record surface value.
logfO2vsIW = zeros(size(r,1), size(r,2));
logfO2 = zeros(size(r,1), size(r,2));
for j = 1:size(r,1)
    for i = 1:length(r)
        [temp_dIW,temp_fO2] = getfO2_H22(P,T,r(j,i),compSheet);
        logfO2vsIW(j,i) = temp_dIW(1);
        logfO2(j,i) = temp_fO2(1);
    end
end


writematrix(t, file, 'Sheet', sheet)
writematrix(r, file, 'WriteMode', 'append', 'Sheet', sheet)

writematrix(t, file, 'Sheet', sheetIW)
writematrix(logfO2vsIW, file, 'WriteMode', 'append', 'Sheet', sheetIW)

writematrix(t, file, 'Sheet', sheetfO2)
writematrix(logfO2, file, 'WriteMode', 'append', 'Sheet', sheetfO2)

figure(1);
subplot(1,2,1);
hold on
box on
plot(t, r_m_Dt(1:100,:), "LineWidth", 1.5)
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

