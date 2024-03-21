% EL
% Feb 2023
% Updated March 2024

% Figure of P vs T and P vs FeRatio
% geotherm data as given by Miyazaki & Korenaga 2019 
% r_eq values calculated from calcFeRatio.m


clear;

Tp = [2500 3000 3500 4000 4500];
P = 0:0.5e9:120e9;
compSheet_earth = 'H04_E';          %sheet in MoleWeights.xlsx to use for composition, need end of accretion value


% READ DATA SHEETS
PV_data = readmatrix('/db/PVcalc.xlsx');
Adiabat_data = readmatrix('\db\geotherms_combo.xlsx');
CompEarth_data = readmatrix('\db\Compositions.xlsx', 'Sheet', compSheet_earth);
MolW_byM_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11', 'Range', 'B2:B13');
MolW_data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11', 'Range', 'C2:C13');

CompEarth = CompEarth_data(:,end);

Tad = zeros(length(P), length(Tp));
PV = zeros(length(P), length(Tp));
r_eq = zeros(length(P), length(Tp));

for i = 1:length(Tp)

    Tad(:,i) = getMOAdiabat(Tp(i), P, Adiabat_data);

    % DETERMINE PV INTEGRATION AS INT(PV)/RT
    PV(:,i) = calcPV(Tad(:,i),P,PV_data);

    % CALCULATE Fe3+/sumFe EQUILIBRIUM RATIO AS FUNCTION OF P,T,dIW
    [r_eq(:,i),dIW] = calcFeRatio(Tad(:,i), P, PV(:,i), CompEarth, MolW_data, MolW_byM_data);
    
end
    
data_adiab = readmatrix('\db\geotherms_combo.xlsx', 'Sheet', 'Sheet1');
P_data = data_adiab(:,1);
Tsol = data_adiab(:,3);
Tliq = data_adiab(:,4);
Tphi = data_adiab(:,5);

%range for post-Cr oxidation = modern day mantle FeO*
r_low_f = 0.02;
r_high_f = 0.06;

%range for pre-Cr oxidation with 8.1% FeO*
r_low_0 = r_low_f + 0.35/8.2;
r_high_0 = r_high_f + 0.35/8.2;
    
figure('Position', [30 30 1000 500]);
subplot(1,2,1)
hold on
box on
colororder('default')
p1 = plot(Tad(:,1), P/1e9, "LineWidth", 1.5);
p2 = plot(Tad(:,2), P/1e9, "LineWidth", 1.5);
p3 = plot(Tad(:,3), P/1e9, "LineWidth", 1.5);
p4 = plot(Tad(:,4), P/1e9, "LineWidth", 1.5);
p5 = plot(Tad(:,5), P/1e9, "LineWidth", 1.5);
p6 = plot(Tsol, P_data, 'k--');
p7 = plot(Tliq, P_data, 'k-');
p8 = plot(Tphi, P_data, 'k:');
ylabel('Pressure (GPa)')
xlabel('Temperature (K)')
ylim([0 120])
yticks([0, 20, 40, 60, 80, 100, 120])
yticklabels(["0","20","40","60","80","100","120","135"])
ax = gca;
ax.YDir = 'reverse';
ax.YColor = [0.15 0.15 0.15];

text(5500, 5, 'a', 'FontSize', 20, 'FontWeight', 'bold')
legend([p1 p2 p3 p4 p5 p6 p7 p8], "2500 K", "3000 K", "3500 K", "4000 K", "4500 K", "solidus", "liquidus", "\phi = 0.4", 'Location', 'southwest')

subplot(1,2,2)
hold on
box on
colororder('default')
x = [r_low_f r_high_f];       %shaded region
y = [135 135];
area(x,y, 'FaceColor', [0.87 0.87 0.87])
x_0 = [r_low_0 r_high_0];
area(x_0,y, 'FaceColor', [0.75 0.75 0.75])
% plot(r_eq(:,1,1), P/1e9, ':', 'Color', "#0072BD", "LineWidth", 1);
% plot(r_eq(:,2,1), P/1e9, ':', "Color", "#D95319", "LineWidth", 1);
% plot(r_eq(:,3,1), P/1e9, ':', "Color", "#EDB120", "LineWidth", 1);
% plot(r_eq(:,4,1), P/1e9, ':', "Color", "#7E2F8E", "LineWidth", 1);
% plot(r_eq(:,5,1), P/1e9, ':', "Color", "#77AC30", "LineWidth", 1);
% plot(r_eq(:,1,2), P/1e9, '--', 'Color', "#0072BD", "LineWidth", 1);
% plot(r_eq(:,2,2), P/1e9, '--', "Color", "#D95319", "LineWidth", 1);
% plot(r_eq(:,3,2), P/1e9, '--', "Color", "#EDB120", "LineWidth", 1);
% plot(r_eq(:,4,2), P/1e9, '--', "Color", "#7E2F8E", "LineWidth", 1);
% plot(r_eq(:,5,2), P/1e9, '--', "Color", "#77AC30", "LineWidth", 1);
plot(r_eq(:,1), P/1e9, 'Color', "#0072BD", "LineWidth", 1.5);
plot(r_eq(:,2), P/1e9, "Color", "#D95319", "LineWidth", 1.5);
plot(r_eq(:,3), P/1e9, "Color", "#EDB120", "LineWidth", 1.5);
plot(r_eq(:,4), P/1e9, "Color", "#7E2F8E", "LineWidth", 1.5);
plot(r_eq(:,5), P/1e9, "Color", "#77AC30", "LineWidth", 1.5);

xlabel('Fe^{3+}/\SigmaFe Ratio')
xlim([0 0.3])
xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
ylim([0 120])
yticks([0, 20, 40, 60, 80, 100, 120])
yticklabels(["","","","","","","",""])
ax = gca;
ax.YDir = 'reverse';
ax.YColor = [0.15 0.15 0.15];
ax.Layer = 'top';
text(0.27, 5, 'b', 'FontSize', 20, 'FontWeight', 'bold')

hold off


% figure(2);              % simplified Fe Ratio graph for AGU poster
% hold on
% box on
% colororder('default')
% % x = [r_low_f r_high_f];       %shaded region
% % y = [135 135];
% % area(x,y, 'FaceColor', [0.87 0.87 0.87])
% % x_0 = [r_low_0 r_high_0];
% % area(x_0,y, 'FaceColor', [0.75 0.75 0.75])
% % plot(r_eq(:,1,1), P/1e9, ':', 'Color', "#0072BD", "LineWidth", 1);
% % plot(r_eq(:,2,1), P/1e9, ':', "Color", "#D95319", "LineWidth", 1);
% % plot(r_eq(:,3,1), P/1e9, ':', "Color", "#EDB120", "LineWidth", 1);
% % plot(r_eq(:,4,1), P/1e9, ':', "Color", "#7E2F8E", "LineWidth", 1);
% % plot(r_eq(:,5,1), P/1e9, ':', "Color", "#77AC30", "LineWidth", 1);
% % plot(r_eq(:,1,2), P/1e9, '--', 'Color', "#0072BD", "LineWidth", 1);
% % plot(r_eq(:,2,2), P/1e9, '--', "Color", "#D95319", "LineWidth", 1);
% % plot(r_eq(:,3,2), P/1e9, '--', "Color", "#EDB120", "LineWidth", 1);
% % plot(r_eq(:,4,2), P/1e9, '--', "Color", "#7E2F8E", "LineWidth", 1);
% % plot(r_eq(:,5,2), P/1e9, '--', "Color", "#77AC30", "LineWidth", 1);
% r1 = plot(r_eq(:,1,3), P/1e9, 'Color', "#0072BD", "LineWidth", 1.5);
% r2 = plot(r_eq(:,2,3), P/1e9, "Color", "#D95319", "LineWidth", 1.5);
% r3 = plot(r_eq(:,3,3), P/1e9, "Color", "#EDB120", "LineWidth", 1.5);
% r4 = plot(r_eq(:,4,3), P/1e9, "Color", "#7E2F8E", "LineWidth", 1.5);
% r5 = plot(r_eq(:,5,3), P/1e9, "Color", "#77AC30", "LineWidth", 1.5);
% 
% legend([r1 r2 r3 r4 r5], "2500 K", "3000 K", "3500 K", "4000 K", "4500 K", 'Location', 'northeast')
% 
% xlabel('Fe^{3+}/\SigmaFe Ratio')
% ylabel('Pressure (GPa)')
% xlim([0 0.3])
% xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
% ylim([0 120])
% yticks([0, 20, 40, 60, 80, 100, 120])
% yticklabels(["0","20","40","60","80","100","120","135"])
% ax = gca;
% ax.YDir = 'reverse';
% ax.YColor = [0.15 0.15 0.15];
% ax.Layer = 'top';
% 
% hold off