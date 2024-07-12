% EL
% August 2022
% Updated July 2024
% 
% Creates figure for metal rain model as Fe ratio vs time, MC results ONLY.
% Uses data from Rain_MC.xlsx, which is created by Ratios_MetalRain_MC.m
%
% Creates Figure 2 in MO redox manuscript.
%
% To output as vector file:
% print(gcf,'-vector','-dsvg',['MC_Tconst','.svg'])

reset = 1;          % change value to 1 for initialization of data.

if reset == 1
    clear;
    
    sheet_r1 = 'H04_const';      %MC results
    model1 = 4;            %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
    
    sheet_r2 = 'N21_const';      %MC results
    model2 = 5;              %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
        
    % the MC results
    data = readmatrix('\db\Rain_MC.xlsx', 'Sheet', sheet_r1);
    r_0p_1 = data(2,:);          %the first row is t
    r_1p_1 = data(3,:);
    r_5p_1 = data(4,:);
    r_25p_1 = data(5,:);
    r_50p_1 = data(6,:);
    r_75p_1 = data(7,:);
    r_95p_1 = data(8,:);
    r_99p_1 = data(9,:);
    r_100p_1 = data(10,:);
        
    % the MC results
    data = readmatrix('\db\Rain_MC.xlsx', 'Sheet', sheet_r2);
    r_0p_2 = data(2,:);          %the first row is t
    r_1p_2 = data(3,:);
    r_5p_2 = data(4,:);
    r_25p_2 = data(5,:);
    r_50p_2 = data(6,:);
    r_75p_2 = data(7,:);
    r_95p_2 = data(8,:);
    r_99p_2 = data(9,:);
    r_100p_2 = data(10,:);

end

[t, Accr_model1] = getAccrModel(model1);
[~, Accr_model2] = getAccrModel(model2);

% 0.35% FeO1.5 destroyed by Cr oxidation, see Hirschmann 2022
r_preCr_MC_1 = [r_0p_1(end), r_1p_1(end), r_5p_1(end), r_25p_1(end), r_50p_1(end), ...
                r_75p_1(end), r_95p_1(end), r_99p_1(end), r_100p_1(end)];
r_postCr_MC_1 = max(r_preCr_MC_1 - 0.35/8.2, 0);

r_preCr_MC_2 = [r_0p_2(end), r_1p_2(end), r_5p_2(end), r_25p_2(end), r_50p_2(end), ...
                r_75p_2(end), r_95p_2(end), r_99p_2(end), r_100p_2(end)];
r_preCr_MC_2_shift = r_preCr_MC_2 + 0.005;          %shifted for clarity in graph
r_postCr_MC_2 = max(r_preCr_MC_2 - 0.35/8.2, 0);

r_preCr_limit = 0.35/8.2;       %zero limit for preCr

r_postCr_range = [0.02 0.06];
r_preCr_range = r_postCr_range + 0.35/8.2;

colors_blue = [33,102,172; 116, 159, 210; 172, 192, 224; 224, 231, 243; 255, 255, 255]/255;
colors_red = [213,62,79; 229, 134, 145; 239, 183, 189; 249, 228, 230; 255, 255, 255]/255;

figure('Units', 'centimeters', 'Position', [3 3 36 13]);

subplot(1,3,1);     %Fe ratio for H04
hold on
box on
t2 = [t, fliplr(t)];
mid50 = [r_25p_1, fliplr(r_75p_1)];
mid90 = [r_5p_1, fliplr(r_95p_1)];
mid99 = [r_1p_1, fliplr(r_99p_1)];
fill(t2, mid99, colors_red(4,:), 'EdgeColor', 'none');
fill(t2, mid90, colors_red(3,:), 'EdgeColor', 'none');
fill(t2, mid50, colors_red(2,:), 'EdgeColor', 'none');
plot(t, r_50p_1, '-', 'Color', colors_red(1,:), "LineWidth", 1.5);
plot(t, r_0p_1, ':', 'Color', colors_red(1,:), "LineWidth", 1);
plot(t, r_100p_1, ':', 'Color', colors_red(1,:), "LineWidth", 1);
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.16])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
xlim([0 100])
text(3, 0.153, 'a', 'FontSize', 20, 'FontWeight', 'bold')
ax = gca;
ax.Layer = 'top';

axes('Position', [0.255 0.175 0.08 0.18]);
box on
plot(t, Accr_model1, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.Layer = 'top';
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
ax.FontSize = 7;
xlabel("Time (Myr)")
ylabel("Mass fraction")


subplot(1,3,2);     %Fe ratio for N21
hold on
box on
t2 = [t, fliplr(t)];
mid50 = [r_25p_2, fliplr(r_75p_2)];
mid90 = [r_5p_2, fliplr(r_95p_2)];
mid99 = [r_1p_2, fliplr(r_99p_2)];
fill(t2, mid99, colors_blue(4,:), 'EdgeColor', 'none');
fill(t2, mid90, colors_blue(3,:), 'EdgeColor', 'none');
fill(t2, mid50, colors_blue(2,:), 'EdgeColor', 'none');
plot(t, r_50p_2, '-', 'Color', colors_blue(1,:), "LineWidth", 1.5);
plot(t, r_0p_2, ':', 'Color', colors_blue(1,:), "LineWidth", 1)
plot(t, r_100p_2, ':', 'Color', colors_blue(1,:), "LineWidth", 1)
%plot(t, r_1p_2, '--', 'Color', colors_blue(1,:))
%plot(t, r_99p_2, '--', 'Color', colors_blue(1,:))
%yline(r_eq(end), '--')
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.16])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16])
xlim([0 100])
text(3, 0.153, 'b', 'FontSize', 20, 'FontWeight', 'bold')
ax = gca;
ax.Layer = 'top';

axes('Position', [0.535 0.175 0.08 0.18]);
box on
plot(t, Accr_model2, 'k', "LineWidth", 1.5)
xticks([0 20 40 60 80 100])
yticks([0 0.2 0.4 0.6 0.8 1.0])
ax = gca;
ax.TickLength = [0.02 0];
ax.XTickLabelRotation = 0;
ax.Layer = 'top';
ax.FontSize = 7;
xlabel("Time (Myr)")
ylabel("Mass fraction")


ax1 = subplot(1,3,3);
hold on
box on

ax2 = axes;
hold on
box on
x_lim = [0.0 0.18];
y_lim = [0 0.12];
x2 = [0, 0, x_lim(2), x_lim(2)];
y2 = [r_postCr_range(1), r_postCr_range(2), r_postCr_range(2), r_postCr_range(1)];
patch('Xdata', x2, 'Ydata', y2, 'FaceColor', [0.87 0.87 0.87], 'EdgeColor', 'none')
c = [5, 4, 3, 2, 1, 2, 3, 4, 5];
%r_preCr_MC_1_wr0 = sort(cat(2, r_preCr_MC_1, r_preCr_limit));           %where r = 0, flatten to line
%r_postCr_MC_1_wr0 = sort(cat(2, r_postCr_MC_1, 0));
plot(r_preCr_MC_1, r_postCr_MC_1, '-', 'Color', colors_red(1,:), "LineWidth", 2)
scatter(r_preCr_MC_1, r_postCr_MC_1, [], c, 'filled', 'MarkerEdgeColor', colors_red(1,:))
xlabel('Fe^{3+}/\SigmaFe pre-Cr oxidation')
ylabel('Fe^{3+}/\SigmaFe post-Cr oxidation')
xlim(x_lim)
ylim(y_lim)
xticks([0 0.03 0.06 0.09 0.12 0.15 0.18])
yticks([0 0.02 0.04 0.06 0.08 0.10 0.12])
ax2.Layer = 'top';
% labels1 = ["P_0   "];
% r_preCr_MC_1_rightshift = r_preCr_MC_1+0.002;
% r_postCr_MC_1_upshift = r_postCr_MC_2+0.004;
% text(r_preCr_MC_1_rightshift(1), r_postCr_MC_1_upshift(1), labels1, 'FontSize', 6, 'Color', colors_red(1,:), 'HorizontalAlignment', 'left');
labels2 = ["P_0   ", "P_1   ", "P_5   ", "P_{25}   ", "median   ", "P_{75}   ", "P_{95}   ", "P_{99}   ", "P_{100}   "];
text(r_preCr_MC_1, r_postCr_MC_1, labels2, 'FontSize', 8, 'Color', colors_red(1,:), 'HorizontalAlignment', 'right');

ax3 = axes;
hold on
box on
% r_preCr_MC_2_wr0 = sort(cat(2, r_preCr_MC_2_shift, r_preCr_limit+0.005));           %where r = 0, flatten to line
r_postCr_MC_2_wr0 = sort(cat(2, r_postCr_MC_2, 0));
plot(r_preCr_MC_2_shift, r_postCr_MC_2, '-', 'Color', colors_blue(1,:), "LineWidth", 2)
scatter(r_preCr_MC_2_shift, r_postCr_MC_2, [], c, 'filled', 'MarkerEdgeColor', colors_blue(1,:))
% labels1 = ["P_0", ];
% r_postCr_MC_2_upshift = r_postCr_MC_2+0.004;
% r_preCr_MC_2_leftshift = r_preCr_MC_2_shift-0.003;
% text(r_preCr_MC_2_leftshift(1), r_postCr_MC_2_upshift(1), labels1, 'FontSize', 6, 'Color', colors_blue(1,:), 'HorizontalAlignment', 'left');
labels2 = ["   P_0", "   P_1", "   P_5", "   P_{25}", "   median", "   P_{75}", "   P_{95}", "   P_{99}", "   P_{100}"];
text(r_preCr_MC_2_shift, r_postCr_MC_2, labels2, 'FontSize', 8, 'Color', colors_blue(1,:), 'HorizontalAlignment', 'left');

text(0.005, 0.115, 'c', 'FontSize', 20, 'FontWeight', 'bold')

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

