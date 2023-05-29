% EL
% August 2022
% Updated 2023-03-12
%
% Determine Fe3/Fe in MO of given depth (d_mo) and given efficiency (eff)
% for late accretion (1%)


clear;

data = readmatrix('geotherms.xlsx', 'Sheet', '2500K');
P = data(:,2)*1e9;
z = data(:,3)*1e3;
T = data(:,4);
r_eq = data(:,7);           %6 = avg, 7 = D20, 8 = A19

%Fe3/sumFe value AFTER
%r_0 = 0.1267;     %median from H04 modeling
r_0 = 0.1116;     %median from N21 modeling
%r_0 = 0.02 + 0.35/8.05;                   %low modern value
%r_0 = 0.06 + 0.35/8.05;                   %high modern value

sheet_nomix = 'N21_nomix';     %sheet name to record data
sheet_mix = 'N21_mix';
file = 'Rain_EffvsD_late.xlsx';               % file name


% CONSTANTS
Mm = 4e24;        %[kg] mass of whole mantle
R_E_0 = 6371e3;   %[m] radius of Earth
rho_E = 5500;     %[kg/m^3] mean density of Earth
rho_Si = 3750;    %[kg/m^3] silicate density
D = 1e-7;         %[m^2/s] chemical diffusity
d = 1000;         %[m] eq distance from Rubie+ 2003 (100 to 1000 m)

M_late = 0.01;      % mass accreted in late accretion
M_prelate = 0.99;   % mass post-GI
M_postlate = M_late+M_prelate;

% estimate time evolution of iron on Earth
M_earth = 5.97e24;                  % [kg] mass of Earth
Fe_earth = 0.321;                   % [] weight % of iron
M_Fe_earth = M_earth*Fe_earth;      % [kg] estimated mass of iron on Earth
M_Fe_accr = M_Fe_earth*M_late;        % [kg] mass of iron in GI

% radius of Earth over time
R_E_prelate = (3*M_prelate*M_earth/(4*pi*rho_E)).^(1/3);
R_E_postlate = (3*(M_prelate+M_late)*M_earth/(4*pi*rho_E)).^(1/3);

% pressure at base of MO (assume base of MO = proportional to
% present day mantle depth)
z_base_prelate = z(end)*R_E_prelate/R_E_0;
z_base_postlate = z(end)*R_E_postlate/R_E_0;
P_base_prelate = interp1(z,P,z_base_prelate);
P_base_postlate = interp1(z,P,z_base_postlate);

% iron droplets
r_d = 0.5e-2;                   %[m] radius of droplets (~1cm diameter)
v_d = 0.5;                      %[m/s] percolation velocity
rho_Fe = 7800;                  %[kg/m^3] iron melt density
m_d = rho_Fe * 4/3*pi*r_d^3;    %[kg] mass of a single droplet
dM_Fe = M_Fe_accr;              %[kg] iron mass accreted per time interval
N_d = dM_Fe/m_d;                % # of droplets per time interval

J = 2*pi*rho_Si*D*r_d;  %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
dt = d/v_d;             %[s] time interval to fall equilibrium distance d
dMp_temp = J*dt*N_d;

% DETERMINE POST-GI FE RATIO
eff = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]/100;
[val,idx] = min(abs(z-z_base_prelate));
z_late = linspace(z(2),round(z_base_prelate,-4), idx-1);

r_m = zeros(1, length(r_eq));
effvsd_nomix = zeros(length(z_late), length(eff));
effvsd_mix = zeros(length(z_late), length(eff));

for k = 1:length(eff)          % for each efficiency

    for j = 1:length(z_late)     % for each depth

        r_m(1) = r_0;          % initiate initial r_m as present mantle ratio

        Mmo_imp = rho_Si * 4/3*pi*(R_E_postlate^3 - (R_E_postlate-z_late(j))^3);      %mass of spherical shell MO
        dMp = min(dMp_temp*eff(k), Mmo_imp);                               %silicate mass eq'd with given efficiency
        
        for i = 2:length(r_eq)        % droplets falling through mantle time
            if z(i) <= z_late(j)
                % ratio changes as the droplets fall through mantle layers until base of MO
                r_m(i) = (dMp*r_eq(i) + (Mmo_imp-dMp)*r_m(i-1))/(Mmo_imp);
            else
                r_m(i) = r_m(i-1);
            end
        end
        
        % final r from mixing redox'd MO with whole mantle
        r_m_mix = (Mmo_imp*r_m(end) + (Mm*M_postlate-Mmo_imp)*r_0)/(Mm*M_postlate);

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


writematrix(effvsd_nomix, file, 'Sheet', sheet_nomix)
writematrix(effvsd_mix, file, 'Sheet', sheet_mix)



% THIS MAY OR MAY NOT WORK :)

c = linspace(0.005, 0.1, 20);
ct = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10];
c0 = [0.02 0.03];
% [C0, h0] = contour(z_late/1e3, log10(eff), effvsd_mix', c0);
% Ni = C0(2,1);
% eff0 = zeros(1,Ni)-3;
% z0 = linspace(0,max(z_late)/1e3,Ni);

figure('Position', [200 200 1000 400]);
subplot(1,2,1)
hold on
box on
[C,h] = contour(z_late/1e3, log10(eff), effvsd_mix'-0.35/8.05, c, 'LineWidth', 2);
clabel(C,h,ct, 'LabelSpacing', 200)
xlabel('Depth (km)')
ylabel('log(Efficiency)')
xlim([0 1000])
title("Whole Mantle (with post-mixing)")

subplot(1,2,2)
hold on
box on
% r0 = [eff0, fliplr(C0(2,2:Ni+1))];
% z0b = [z0, fliplr(z0)];
% fill(z0b, r0, [220/255 220/255 220/255]);
%plot(C0(1,2:Ni+1), C0(2,2:Ni+1), 'Color', [200 200 200])
colormap cool
colormap(flipud(colormap));
contour(z_late/1e3, log10(eff), effvsd_nomix'-0.35/8.05, c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 300)
xlabel('Depth (km)')
ylabel('log(Efficiency)')
xlim([0 1000])
ylim([-3 0])
title("Top MO only (no post-mixing)")


