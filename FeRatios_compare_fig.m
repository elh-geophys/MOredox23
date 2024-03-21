% EL
% April 2023
% Updated March 5, 2024

% compares Fe3+/sumFe ratio for different compositions


clear;

dIW = -2;

data = readmatrix('\db\geotherms_combo.xlsx');
P = data(:,1)*1e9;
z = data(:,2);
Tsol = data(:,3);
Tliq = data(:,4);
Tphi = data(:,5);
T = data(:,7:11);

PV_data = readmatrix('/db/PVcalc.xlsx');

comp_D20 = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Deng20');
comp_A19 = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Armstrong19');

[val, cutoff_2500K] = min(abs(Tphi-T(:,1)));

r_eq_D20 = zeros(size(T));
r_eq_A19 = zeros(size(T));
PV_term = zeros(size(T));

% DETERMINE PV INTEGRATION AS INT(PV)/RT
for i = 1:size(T,2)
    PV_term(:,i) = calcPV(T(:,i),P,PV_data);
end

% DETERMINE Fe3+/sumFe equilibrium ratio
for i = 1:size(T,2)
    r_eq_D20(:,i) = calcFeRatio_byIW(T(:,i), P, dIW, PV_term(:,i), comp_D20(:,2), comp_D20(:,3));
    r_eq_A19(:,i) = calcFeRatio_byIW(T(:,i), P, dIW, PV_term(:,i), comp_A19(:,2), comp_A19(:,3));
end

% to get correct scale for depth axis
z_test = [0 500 1000 1500 2000 2500 2800];
P_test = zeros(1, length(z_test));
T_test = zeros(1, length(z_test));
r_test = zeros(1, length(z_test));
for i = 1:length(z_test)
    [~,idx] = min(abs(z - z_test(i)));
    P_test(i) = P(idx);
    T_test(i) = Tphi(idx);
    r_test(i) = r_eq_D20(i,1);
end

P = P/1e9;
    
figure('Position', [30 30 400 500]);
hold on
box on
yyaxis left
colororder('default')
r1 = plot(r_eq_D20(1:cutoff_2500K,1), P(1:cutoff_2500K), 'Color', "#0072BD", "LineWidth", 1.5);
r2 = plot(r_eq_D20(:,2), P, "Color", "#D95319", "LineWidth", 1.5);
r3 = plot(r_eq_D20(:,3), P, "Color", "#EDB120", "LineWidth", 1.5);
r4 = plot(r_eq_D20(:,4), P, "Color", "#7E2F8E", "LineWidth", 1.5);
r5 = plot(r_eq_D20(:,5), P, "-", "Color", "#77AC30", "LineWidth", 1.5);
plot(r_eq_A19(1:cutoff_2500K,1), P(1:cutoff_2500K), ":", 'Color', "#0072BD", "LineWidth", 1.5);
plot(r_eq_A19(:,2), P, ":", "Color", "#D95319", "LineWidth", 1.5);
plot(r_eq_A19(:,3), P, ":", "Color", "#EDB120", "LineWidth", 1.5);
plot(r_eq_A19(:,4), P, ":", "Color", "#7E2F8E", "LineWidth", 1.5);
plot(r_eq_A19(:,5), P, ":", "Color", "#77AC30", "LineWidth", 1.5);
xlabel('Fe^{3+}/\SigmaFe Ratio')
ylabel('Pressure (GPa)')

xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
ylim([0 135])
yticks([0, 20, 40, 60, 80, 100, 120, 135])
yticklabels(["0","20","40","60","80","100","120","135"])
ax = gca;
ax.YDir = 'reverse';
ax.YColor = [0.15 0.15 0.15];
ax.Layer = 'top';

yyaxis right
plot(r_test, P_test/1e9, 'k.', 'MarkerSize', 0.1)
ylabel("Depth (km)")
ylim([0 135.14])
yticks(P_test/1e9)
yticklabels(z_test)
ax = gca;
ax.YDir = 'reverse';
ax.YColor = [0.15 0.15 0.15];

legend([r1 r2 r3 r4 r5], "2500 K", "3000 K", "3500 K", "4000 K", "4500 K", 'Location', 'southwest')


hold off

