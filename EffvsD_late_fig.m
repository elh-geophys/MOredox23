% EL
% Aug 2022
% Updated 2023-03-12
%
% Efficiency vs Depth figure (Late Accretion)

reset = 1;

if reset == 1
    clear;
    
    data = readmatrix('\db\Rain_EffvsD_late.xlsx', 'Sheet', 'data');
    Z = data(1:100,1);
    eff = data(10:28,2);            %just 1%-100%

    %these are all set before Cr oxi    
    data = readmatrix('Rain_EffvsD_late.xlsx', 'Sheet', 'H04_nomix');
    H04_nomix = data(1:100,10:28);
    data = readmatrix('Rain_EffvsD_late.xlsx', 'Sheet', 'H04_mix');
    H04_mix = data(1:100,10:28);
    
    data = readmatrix('Rain_EffvsD_late.xlsx', 'Sheet', 'N21_nomix');
    N21_nomix = data(1:100,10:28);
    data = readmatrix('Rain_EffvsD_late.xlsx', 'Sheet', 'N21_mix');
    N21_mix = data(1:100,10:28);
    
end

mix = H04_mix - 0.35/8.05;       %0.35% reduction in FeO1.5 after Cr oxi with 8.05% FeO* from Deng20 comp
nomix = H04_nomix - 0.35/8.05;
val = 1;                            %for r_0 to choose
title_name = "H04";

%Fe3/sumFe value AFTER GI
r_0_H04 = 0.1267 - 0.35/8.05;     %median from H04 modeling
r_0_N21 = 0.1116 - 0.35/8.05;     %median from N21 modeling

r_0 = [r_0_H04, r_0_N21];
     
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
%map = [0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1; 1 1 1; 1 1 1];
map = [0.72 0.72 0.72; 1 1 1; 1 1 1; 1 1 1; 1 1 1];
%map = [0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];                                       % chosen mix contour colors

map_neg = [0 0 0];

c = linspace(0.005, 0.1, 20);
ct = [0.01 0.02 0.03 0.04 0.05 0.06];
c0 = [0.02 0.03 0.04 0.05 0.06];
c_neg = [-0.03 0];

figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, Z/1e3, log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, Z/1e3, log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');

ax3 = axes;
contour(ax3, Z/1e3, log10(eff), nomix', c_neg, 'LineWidth', 1, 'LineStyle', '--', 'FaceColor', 'none');

linkaxes([ax1,ax2,ax3])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];

colormap(ax1, map)
colormap(ax2, flipud(colormap(ax2,cool)));
colormap(ax3, map_neg);

ax1.Box = 'on';
ax1.Layer = 'top';
ax1.XLabel.String = 'Depth (km)';
ax1.XTick = [0 200 400 600 800 1000];
ax1.YLabel.String = 'Degree of Chemical Equilibration';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};
text(820, -1.9, "r_0=" + round(r_0(val),4), 'FontWeight', 'bold')
title(ax1, title_name)




