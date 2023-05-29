% EL
% August 2022
% Updated 03-07-2023
%
% Determine Fe3/sumFe in MO as a function of given depth (d_mo) and given
% efficiency (eff) for the last giant impact (GI)

clear;

data = readmatrix('\db\geotherms.xlsx', 'Sheet', '3500K');
P = data(:,2)*1e9;
z = data(:,3)*1e3;
T = data(:,4);
r_eq = data(:,7);           %6 = avg, 7 = D20, 8 = A19

%Fe3/sumFe value BEFORE GI from H04 modeling
% r0 = 0.0931;        % 0th
% r0 = 0.1040;        % 1st
% r0 = 0.1101;        % 5th      near upper limit
% r0 = 0.1197;        % 25th
% r0 = 0.1267;        % 50th
% r0 = 0.1338;        % 75th
% r0 = 0.1408;        % 95th
% r0 = 0.1450;        % 99th
% r0 = 0.1500;        % 100th

%Fe3/sumFe value BEFORE GI from N21 modeling
% r0 = 0.0067;        % 0th      impossible
% r0 = 0.0489;        % 1st      works
% r0 = 0.0738;        % 5th
% r0 = 0.0990;        % 25th
% r0 = 0.1116;        % 50th     above upper limit (0.1102)
% r0 = 0.1220;        % 75th
% r0 = 0.1321;        % 95th
% r0 = 0.1357;        % 99th
% r0 = 0.1402;        % 100th

%test value
%r_0 = 0.06 + 0.35/8.05;
%r_0 = 0.02 + 0.35/8.05;
r_0 = 0.095;

sheet_nomix = 'H04_25th_nomix';     %sheet name to record data
sheet_mix = 'H04_25th_mix';
file = 'Rain_EffvsD_MConly.xlsx';        % file name
write = 0;                          %1 to write, else to not write to file

% CHANGE GI VALUE TOO!
% H04
%M_GI = 0.09;      % mass accreted in GI
%M_preGI = 0.90;   % mass pre-GI
% N21
M_GI = 0.49;      % mass accreted in GI
M_preGI = 0.50;   % mass pre-GI

M_postGI = M_GI+M_preGI;


% CONSTANTS
Mm = 4e24;        %[kg] mass of whole mantle
R_E_0 = 6371e3;   %[m] radius of Earth
rho_E = 5500;     %[kg/m^3] mean density of Earth
rho_Si = 3750;    %[kg/m^3] silicate density
D = 1e-7;         %[m^2/s] chemical diffusity
d = 1000;         %[m] eq distance from Rubie+ 2003 (100 to 1000 m)

% estimate time evolution of iron on Earth
M_earth = 5.97e24;                  % [kg] mass of Earth
Fe_earth = 0.321;                   % [] weight % of iron
M_Fe_earth = M_earth*Fe_earth;      % [kg] estimated mass of iron on Earth
M_Fe_accr = M_Fe_earth*M_GI;        % [kg] mass of iron in GI

% radius of Earth over time
R_E_preGI = (3*M_preGI*M_earth/(4*pi*rho_E)).^(1/3);
R_E_postGI = (3*(M_preGI+M_GI)*M_earth/(4*pi*rho_E)).^(1/3);

% pressure at base of MO (assume base of MO = proportional to
% present day mantle depth)
z_base_preGI = z(end)*R_E_preGI/R_E_0;
z_base_postGI = z(end)*R_E_postGI/R_E_0;
P_base_preGI = interp1(z,P,z_base_preGI);
P_base_postGI = interp1(z,P,z_base_postGI);

% iron droplets
r_d = 0.5e-2;                   %[m] radius of droplets (~1cm diameter)
v_d = 0.5;                      %[m/s] percolation velocity
rho_Fe = 7800;                  %[kg/m^3] iron melt density
m_d = rho_Fe * 4/3*pi*r_d^3;    %[kg] mass of a single droplet
dM_Fe = M_Fe_accr;    %[kg] iron mass accreted per time interval
N_d = dM_Fe/m_d;                % # of droplets per time interval

J = 2*pi*rho_Si*D*r_d;  %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
dt = d/v_d;             %[s] time interval to fall equilibrium distance d
dMp_temp = J*dt*N_d;

% DETERMINE POST-GI FE RATIO
eff = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]/100;
[val,idx] = min(abs(z-z_base_preGI));
z_GI = linspace(z(2),round(z_base_preGI,-4), idx-1);

r_m = zeros(1, length(r_eq));
effvsd_nomix = zeros(length(z_GI), length(eff));
effvsd_mix = zeros(length(z_GI), length(eff));

for k = 1:length(eff)          % for each efficiency

    for j = 1:length(z_GI)     % for each depth

        r_m(1) = r_0;          % initiate initial r_m as present mantle ratio

        Mmo_imp = rho_Si * 4/3*pi*(R_E_postGI^3 - (R_E_postGI-z_GI(j))^3);      %mass of spherical shell MO
        dMp = min(dMp_temp*eff(k), Mmo_imp);                               %silicate mass eq'd with given efficiency
        
        for i = 2:length(r_eq)        % droplets falling through mantle time
            if z(i) <= z_GI(j)
                % ratio changes as the droplets fall through mantle layers until base of MO
                r_m(i) = (dMp*r_eq(i) + (Mmo_imp-dMp)*r_m(i-1))/(Mmo_imp);
            else
                r_m(i) = r_m(i-1);
            end
        end
        
        % final r from mixing redox'd MO with whole mantle
        r_m_mix = (Mmo_imp*r_m(end) + (Mm*M_postGI-Mmo_imp)*r_0)/(Mm*M_postGI);

        effvsd_nomix(j,k) = r_m(end);
        effvsd_mix(j,k) = r_m_mix;
    end

end

% determine logfO2 as dIW for the MO, record surface value.
% logfO2vsIW = zeros(1, length(r_m_Dt));
% for i = 1:length(r_m_Dt)
%     temp = fO2_EL(P/1e9, T, dV, r_m_Dt(i), G_flag, PV_method);
%     logfO2vsIW(i) = temp(1);
% end

if write == 1
    writematrix(effvsd_nomix, file, 'Sheet', sheet_nomix)
    writematrix(effvsd_mix, file, 'Sheet', sheet_mix)
end

%range for post-Cr oxidation = modern day mantle FeO*
r_low_f = 0.02;
r_high_f = 0.06;

%range for pre-Cr oxidation with 8.05% FeO* from Deng20 composition
r_low_0 = r_low_f + 0.35/8.05;
r_high_0 = r_high_f + 0.35/8.05;

%contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16];
c0 = [0.0635 0.0735 0.0835 0.0935 0.1035];
%c0 = [0.07 0.075 0.080 0.085 0.090];           %3-4.5% Fe3/sumFe

figure('Position', [200 200 500 400]);
ax1 = axes;
contourf(ax1, z(1:length(z_GI))/1e3, log10(eff), effvsd_mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, z(1:length(z_GI))/1e3, log10(eff), effvsd_nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

map = [0.88 0.88 0.88; 0.80 0.80 0.80; 0.72 0.72 0.72; 1 1 1; 1 1 1];   %use depending on # of contours
%map = [0.72 0.72 0.72; 1 1 1; 1 1 1];                                      
colormap(ax1, map)
colormap(ax2, flipud(colormap(ax2,cool)));

ax1.Box = 'on';
ax1.XLabel.String = 'Depth (km)';
ax1.YLabel.String = 'Efficiency';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'0.1%' '' '' '' '' '' '' '' '' ...
    '1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};


