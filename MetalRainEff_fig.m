% EL
% August 2022
% Updated March 2023
% 
% Creates figure for metal rain model, as Fe ratio vs time with varying 
% percent values. Also overlays MC results.
% Uses data from Rain_***.xlsx, which is created by Ratios_MetalRain_***.m
%
% Creates Figure 2 in MO_redox manuscript.

reset = 1;          % change value to 1 for initialization of data.

if reset == 1
    clear;
    
    xlsx1 = '\db\Rain_H04.xlsx';
    sheet_eff1 = 'H04';      %fixed efficiencies
    sheet_r1 = 'H04_r';      %MC results
    model1 = 4;              %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
    
    xlsx2 = '\db\Rain_N21.xlsx';
    sheet_eff2 = 'N21';      %fixed efficiencies
    sheet_r2 = 'N21_r';      %MC results
    model2 = 5;              %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
    
    data = readmatrix('\db\geotherms.xlsx', 'Sheet', '3500K');
    r_eq = data(:,7);        %6 = avg, 7 = D20, 8 = A19
    
    %choose 100%
    row_r100 = 2;
    row_fO2100 = 3;
    
    %choose 10%
    row_r10 = 4;
    row_fO210 = 5;
    
    %choose 1%
    row_r1 = 6;
    row_fO21 = 7;
    
    %choose 0.1%
    row_rp1 = 8;
    row_fO2p1 = 9;

    % efficiency fixed values
    data = readmatrix(xlsx1, 'Sheet', sheet_eff1);
    t = data(1,2:end);
    r100_1 = data(row_r100,2:end);
    r10_1 = data(row_r10,2:end);
    r1_1 = data(row_r1,2:end);
    rp1_1 = data(row_rp1,2:end);
    
    data = readmatrix(xlsx2, 'Sheet', sheet_eff2);
    r100_2 = data(row_r100,2:end);
    r10_2 = data(row_r10,2:end);
    r1_2 = data(row_r1,2:end);
    rp1_2 = data(row_rp1,2:end);
    
    % the MC results
    data = readmatrix('\db\Rain_MC.xlsx', 'Sheet', sheet_r1);
    r_avg_1 = data(2,2:end);          %the first row is t
    r_5p_1 =  data(3,2:end);
    r_25p_1 = data(4,2:end);
    r_75p_1 = data(5,2:end);
    r_95p_1 = data(6,2:end);
    
    % the MC results
    data = readmatrix('\db\Rain_MC.xlsx', 'Sheet', sheet_r2);
    r_avg_2 = data(2,2:end);          %the first row is t
    r_5p_2 =  data(3,2:end);
    r_25p_2 = data(4,2:end);
    r_75p_2 = data(5,2:end);
    r_95p_2 = data(6,2:end);

end

[t, Accr_model1] = getAccrModel(model1);
[~, Accr_model2] = getAccrModel(model2);

% 0.35% FeO1.5 destroyed by Cr oxidation, see Hirschmann 2022
r_preCr_eff_1 = [r100_1(end), r10_1(end), r1_1(end), rp1_1(end)];
r_preCr_MC_1 = [r_avg_1(end), r_5p_1(end), r_25p_1(end), r_75p_1(end), r_95p_1(end)];

r_postCr_eff_1 = r_preCr_eff_1 - 0.35/8.05;
r_postCr_MC_1 = r_preCr_MC_1 - 0.35/8.05;

r_preCr_eff_2 = [r100_2(end), r10_2(end), r1_2(end), rp1_2(end)];
r_preCr_MC_2 = [r_avg_2(end), r_5p_2(end), r_25p_2(end), r_75p_2(end), r_95p_2(end)];
r_preCr_MC_2_shift = r_preCr_MC_2 + 0.005;

r_postCr_eff_2 = r_preCr_eff_2 - 0.35/8.05;
r_postCr_MC_2 = r_preCr_MC_2 - 0.35/8.05;

r_postCr_range = [0.02 0.06];
r_preCr_range = r_postCr_range + 0.35/8.05;

colors = 1/255*[0 0 0; ...
                80 80 80; ...
                120 120 120; ...
                160 160 160; ...
                200 200 200; ...
                210 210 210];

colors_blue = [33,102,172; 133, 171, 210; 214, 227, 239]/255;
colors_red = [213,62,79; 235, 164, 173; 248, 219, 222]/255;

figure('Position', [30 30 1400 500]);

subplot(1,3,1);     %Fe ratio for H04
hold on
box on
t2 = [t, fliplr(t)];
mid50 = [r_25p_1, fliplr(r_75p_1)];
mid90 = [r_5p_1, fliplr(r_95p_1)];
fill(t2, mid90, colors_red(3,:), 'EdgeColor', 'none');
fill(t2, mid50, colors_red(2,:), 'EdgeColor', 'none');
avg = plot(t, r_avg_1, '-', 'Color', colors_red(1,:), "LineWidth", 1.5);
plot(t, r_5p_1, '-', 'Color', colors_red(1,:))
plot(t, r_25p_1, '-', 'Color', colors_red(1,:))
plot(t, r_75p_1, '-', 'Color', colors_red(1,:))
plot(t, r_95p_1, '-', 'Color', colors_red(1,:))
Fe1 = plot(t, r100_1, 'Color', colors(1,:), "LineWidth", 1.5);
Fe2 = plot(t, r10_1, 'Color', colors(2,:), "LineWidth", 1.5);
Fe3 = plot(t, r1_1, 'Color', colors(3,:), "LineWidth", 1.5);
Fe4 = plot(t, rp1_1, 'Color', colors(4,:), "LineWidth", 1.5);
yline(r_eq(end), '--')
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.18])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18])
xlim([0 100])
text(90, 0.173, 'a', 'FontSize', 20, 'FontWeight', 'bold')
legend([Fe1 Fe2 Fe3 Fe4], "100%", "10%", "1%", "0.1%", "Location", "southwest")

axes('Position', [0.25 0.19 0.08 0.18]);
box on
plot(t, Accr_model1, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.Layer = 'top';
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
xlabel("Time (Myr)")
ylabel("Mass fraction")


subplot(1,3,2);     %Fe ratio for N21
hold on
box on
t2 = [t, fliplr(t)];
mid50 = [r_25p_2, fliplr(r_75p_2)];
mid90 = [r_5p_2, fliplr(r_95p_2)];
fill(t2, mid90, colors_blue(3,:), 'EdgeColor', 'none');
fill(t2, mid50, colors_blue(2,:), 'EdgeColor', 'none');
avg = plot(t, r_avg_2, '-', 'Color', colors_blue(1,:), "LineWidth", 1.5);
plot(t, r_5p_2, '-', 'Color', colors_blue(1,:))
plot(t, r_25p_2, '-', 'Color', colors_blue(1,:))
plot(t, r_75p_2, '-', 'Color', colors_blue(1,:))
plot(t, r_95p_2, '-', 'Color', colors_blue(1,:))
Fe1 = plot(t, r100_2, 'Color', colors(1,:), "LineWidth", 1.5);
Fe2 = plot(t, r10_2, 'Color', colors(2,:), "LineWidth", 1.5);
Fe3 = plot(t, r1_2, 'Color', colors(3,:), "LineWidth", 1.5);
Fe4 = plot(t, rp1_2, 'Color', colors(4,:), "LineWidth", 1.5);
yline(r_eq(end), '--')
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.18])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16 0.18])
xlim([0 100])
text(90, 0.173, 'b', 'FontSize', 20, 'FontWeight', 'bold')
%legend([Fe1 Fe2 Fe3 Fe4 avg], "100%", "10%", "1%", "0.1%", "median", "Location", "southwest")

axes('Position', [0.53 0.19 0.08 0.18]);
box on
plot(t, Accr_model2, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
ax.Layer = 'top';
xlabel("Time (Myr)")
ylabel("Mass fraction")


ax1 = subplot(1,3,3);
hold on
box on

ax2 = axes;
hold on
box on
x_lim = [0 0.18];
y_lim = [0 0.10];
x2 = [0, 0, x_lim(2), x_lim(2)];
y2 = [r_postCr_range(1), r_postCr_range(2), r_postCr_range(2), r_postCr_range(1)];
patch('Xdata', x2, 'Ydata', y2, 'FaceColor', [0.87 0.87 0.87], 'EdgeColor', 'none')
c2 = [1, 3, 2, 2, 3];
plot(r_preCr_MC_1, r_postCr_MC_1, '-', 'Color', colors_red(1,:))
scatter(r_preCr_MC_1, r_postCr_MC_1, [], c2, 'filled', 'MarkerEdgeColor', colors_red(1,:))
xlabel('Fe^{3+}/\SigmaFe pre-Cr oxidation')
ylabel('Fe^{3+}/\SigmaFe post-Cr oxidation')
xlim(x_lim)
ylim(y_lim)
xticks([0 0.03 0.06 0.09 0.12 0.15 0.18])
yticks([0 0.02 0.04 0.06 0.08 0.10])
ax2.Layer = 'top';
labels = ["median   ", "P_5   ", "P_{25}   ", "P_{75}   ", "P_{95}   "];
text(r_preCr_MC_1, r_postCr_MC_1, labels, 'FontSize', 8, 'Color', colors_red(1,:), 'HorizontalAlignment', 'right');

ax3 = axes;
hold on
box on
plot(r_preCr_MC_2_shift, r_postCr_MC_2, '-', 'Color', colors_blue(1,:))
scatter(r_preCr_MC_2_shift, r_postCr_MC_2, [], c2, 'filled', 'MarkerEdgeColor', colors_blue(1,:))
labels = ["   median", "    P_5", "   P_{25}", "   P_{75}", "   P_{95}"];
text(r_preCr_MC_2_shift, r_postCr_MC_2, labels, 'FontSize', 8, 'Color', colors_blue(1,:), 'HorizontalAlignment', 'left');
text(0.162, 0.096, 'c', 'FontSize', 20, 'FontWeight', 'bold')

linkaxes([ax2,ax3])
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

linkprop([ax1 ax2 ax3],'Position');

colormap(ax2, colors_red)
colormap(ax3, colors_blue);

