% EL
% Aug 2022
% Updated 2023-10-03
%
% Efficiency vs Depth figure

reset = 0;      %reset=1 will read the file again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_Tconst.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'H04_data', 'Range', 1);
    P = data(2,:)/1e9;
    eff = data(1,1:28);
 
%     H04_p5_nomix = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
%     H04_p5_mix = readmatrix(xlsx, 'Sheet', 'H04_5th_mix'); 
    
    H04_p25_nomix = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_p25_mix = readmatrix(xlsx, 'Sheet', 'H04_25th_mix'); 
    
    H04_p50_nomix = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_p50_mix = readmatrix(xlsx, 'Sheet', 'H04_50th_mix'); 
    
    H04_p75_nomix = readmatrix(xlsx, 'Sheet', 'H04_75th_nomix');
    H04_p75_mix = readmatrix(xlsx, 'Sheet', 'H04_75th_mix'); 
    
%     N21_p5_nomix = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
%     N21_p5_mix = readmatrix(xlsx, 'Sheet', 'N21_5th_mix'); 
    
    N21_p25_nomix = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_p25_mix = readmatrix(xlsx, 'Sheet', 'N21_25th_mix'); 
    
    N21_p50_nomix = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
    N21_p50_mix = readmatrix(xlsx, 'Sheet', 'N21_50th_mix'); 
    
    N21_p75_nomix = readmatrix(xlsx, 'Sheet', 'N21_75th_nomix');
    N21_p75_mix = readmatrix(xlsx, 'Sheet', 'N21_75th_mix');
  
end

mix = H04_p75_mix-0.35/8.1;       %0.35% reduction in FeO1.5 after Cr oxi with 8.1% FeO* from Deng20 comp
nomix = H04_p75_nomix-0.35/8.1;
r0_val = 6;
title_name = "H04: 75th Percentile";
letter = "c";

%Fe3/sumFe value BEFORE GI from modeling
%     [0th    1st    5th    25th   50th   75th   95th   99th   100th]
%Tconst method
    r_0 = [0.0740,0.0786,0.0852,0.0970,0.1050,0.1125,0.1216,0.1244,0.1263];     %H04
    %r_0 = [0.0241,0.0370,0.0422,0.0499,0.0552,0.0596,0.0669,0.0705,0.0719];     %N21
%Pmo method
    %r_0 = [0.0466,0.0593,0.0727,0.0964,0.1123,0.1270,0.1434,0.1461,0.1498];     %H04
    %r_0 = [0.0151,0.0184,0.0214,0.0312,0.0422,0.0544,0.0694,0.0744,0.0761];     %N21
%U2Q method
    %r_0 = [0.0345,0.0392,0.0450,0.0514,0.0565,0.0634,0.0709,0.0752,0.0791];     %H04
    %r_0 = [0.0129,0.0175,0.0203,0.0290,0.0371,0.0479,0.0599,0.0629,0.0648];     %N21

%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
%map = [0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
%map = [0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1; 1 1 1];
map = [0.72 0.72 0.72; 1 1 1; 1 1 1];                                       % chosen mix contour colors

map_neg = [0 0 0];

eff_limit_low = 0.30;               % From Zube+2019, 30% efficiency for Hf-W
eff_limit_high = 0.70;              % Z19, 70% for fast accretion (GT)

%range for post-Cr oxidation = modern day mantle FeO*
%0.02 to 0.06 contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16];
c0 = [0.02 0.03 0.04 0.05 0.06];
c_neg = [-0.05 0];                                  %for negative Fe3/sumFe limit
c_r0 = [-0.1, r_0(r0_val)-0.35/8.1];                %for original Fe3/sumFe value

red = 1/255*[213,62,79];

figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, P(1:size(mix,1)), log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, P(1:size(nomix,1)), log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');
yline(log10(eff_limit_low), '-', 'Color', red, 'LineWidth', 1.5);
yline(log10(eff_limit_high), '--', 'Color', red, 'LineWidth', 1.5);

ax3 = axes;
contour(ax3, P(1:size(nomix,1)), log10(eff), nomix', c_neg, 'LineWidth', 2, 'LineStyle', ':', 'FaceColor', 'none');

ax4 = axes;
contour(ax4, P(1:size(nomix,1)), log10(eff), nomix', c_r0, 'LineWidth', 1, 'LineStyle', '--', 'Color', 'b', 'FaceColor', 'none');

linkaxes([ax1,ax2,ax3,ax4])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
ax4.Visible = 'off';
ax4.XTick = [];
ax4.YTick = [];

colormap(ax1, map)
colormap(ax2, flipud(colormap(ax2,cool)));
colormap(ax3, map_neg);

ax1.Box = 'on';
ax1.Layer = 'top';
ax1.XLabel.String = 'Pressure (GPa)';
ax1.YLabel.String = 'Degree of Chemical Equilibration';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'0.1%' '' '' '' '' '' '' '' '' ...
    '1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};
%text(5, -2.8, "r_0=" + round(r_0(val),3), 'FontWeight', 'bold')
%text(72, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
%text(72, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)
text(108, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
text(118, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)
title(ax1, title_name)
