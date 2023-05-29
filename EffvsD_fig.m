% EL
% Aug 2022
% Updated 2023-03-12
%
% Efficiency vs Depth figure

reset = 1;      %reset=1 will read the file again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_MConly.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'data');
    Z = data(1:279,1);
    eff = data(1:28,2);

    H04_p1_nomix = readmatrix(xlsx, 'Sheet', 'H04_1st_nomix');
    H04_p1_mix = readmatrix(xlsx, 'Sheet', 'H04_1st_mix');
    
    H04_p5_nomix = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    H04_p5_mix = readmatrix(xlsx, 'Sheet', 'H04_5th_mix'); 
    
    H04_p25_nomix = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_p25_mix = readmatrix(xlsx, 'Sheet', 'H04_25th_mix'); 
   
    N21_p1_nomix = readmatrix(xlsx, 'Sheet', 'N21_1st_nomix');
    N21_p1_mix = readmatrix(xlsx, 'Sheet', 'N21_1st_mix');
    
    N21_p5_nomix = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
    N21_p5_mix = readmatrix(xlsx, 'Sheet', 'N21_5th_mix'); 
    
    N21_p25_nomix = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_p25_mix = readmatrix(xlsx, 'Sheet', 'N21_25th_mix'); 
  
end

mix = N21_p5_mix - 0.35/8.05;       %0.35% reduction in FeO1.5 after Cr oxi with 8.05% FeO* from Deng20 comp
nomix = N21_p5_nomix - 0.35/8.05;
val = 2;
title_name = "N21: 5th Percentile";

%Fe3/sumFe value BEFORE GI from modeling, as [low90% low50% median high50% high90%]
%r_0 = [0.1040, 0.1101, 0.1197] - 0.35/8.05;     %H04
r_0 = [0.0489, 0.0738, 0.0990] - 0.35/8.05;     %N21

%Fe3/sumFe value BEFORE GI from H04 modeling
% r0 = 0.1040;        % 1st
% r0 = 0.1101;        % 5th      near upper limit
% r0 = 0.1197;        % 25th

%Fe3/sumFe value BEFORE GI from N21 modeling
% r0 = 0.0489;        % 1st      works
% r0 = 0.0738;        % 5th
% r0 = 0.0990;        % 25th

%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
map = [0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1]; 
%map = [0.72 0.72 0.72; 1 1 1; 1 1 1];                                       % chosen mix contour colors

eff_limit_low = 0.30;               % From Zube+2019, 30% efficiency for Hf-W
eff_limit_high = 0.70;              % Z19, 70% for fast accretion (GT)

%range for post-Cr oxidation = modern day mantle FeO*
r_low_f = 0.02;
r_high_f = 0.06;

%contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16];
c0 = [0.02 0.03 0.04 0.05 0.06];

red = 1/255*[213,62,79];

figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, Z(1:size(mix,1))/1e3, log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, Z(1:size(nomix,1))/1e3, log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');
yline(log10(eff_limit_low), '-', 'Color', red, 'LineWidth', 1.5);
yline(log10(eff_limit_high), '--', 'Color', red, 'LineWidth', 1.5);

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

colormap(ax1, map)
colormap(ax2, flipud(colormap(ax2,cool)));

ax1.Box = 'on';
ax1.Layer = 'top';
ax1.XLabel.String = 'Depth (km)';
ax1.YLabel.String = 'Degree of Chemical Equilibration';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'0.1%' '' '' '' '' '' '' '' '' ...
    '1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};
%text(2350, -2.8, "r_0=" + round(r_0(val),3), 'FontWeight', 'bold')
text(1950, -2.8, "r_0=" + round(r_0(val),3), 'FontWeight', 'bold')
title(ax1, title_name)
