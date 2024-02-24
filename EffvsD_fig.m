% EL
% Aug 2022
% Updated 2024-02-22
%
% Efficiency vs Depth for last GI figure
%
% Note: need to change the data sheet because H04 and N21 have different
% ranges of pressure (line 17)

reset = 0;      %reset=1 will read the file again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_Pmo.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'N21_data', 'Range', 1);
    P = data(2,:)/1e9;
    eff = data(1,1:28);
 
    % H04_p5_nomix = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    % H04_p5_mix = readmatrix(xlsx, 'Sheet', 'H04_5th_mix'); 
    % 
    H04_p25_nomix = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_p25_mix = readmatrix(xlsx, 'Sheet', 'H04_25th_mix'); 
    
    H04_p50_nomix = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_p50_mix = readmatrix(xlsx, 'Sheet', 'H04_50th_mix'); 
    
    H04_p75_nomix = readmatrix(xlsx, 'Sheet', 'H04_75th_nomix');
    H04_p75_mix = readmatrix(xlsx, 'Sheet', 'H04_75th_mix'); 
    % 
    % N21_p5_nomix = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
    % N21_p5_mix = readmatrix(xlsx, 'Sheet', 'N21_5th_mix'); 
    % 
    N21_p25_nomix = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_p25_mix = readmatrix(xlsx, 'Sheet', 'N21_25th_mix'); 

    N21_p50_nomix = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
    N21_p50_mix = readmatrix(xlsx, 'Sheet', 'N21_50th_mix'); 

    N21_p75_nomix = readmatrix(xlsx, 'Sheet', 'N21_75th_nomix');
    N21_p75_mix = readmatrix(xlsx, 'Sheet', 'N21_75th_mix');
    % 
end

mix = N21_p25_mix-0.35/8.2;       %0.35% reduction in FeO1.5 after Cr oxi with 8.2% FeO* from Deng20 comp
nomix = N21_p25_nomix-0.35/8.2;
r0_val = 4;
title_name = "N21: 25th Percentile";
letter = "d";

%Fe3/sumFe value BEFORE GI from modeling
%     [0th    1st    5th    25th   50th   75th   95th   99th   100th]
%Tconst method
    %r_0 = [0.0651,0.0794,0.0867,0.1001,0.1085,0.1170,0.1263,0.1297,0.1314];      %H04
    %r_0 = [0.0317,0.0381,0.0448,0.0533,0.0591,0.0629,0.0690,0.0724,0.0742];      %N21
%Pmo method
    %r_0 = [0.0415,0.0533,0.0680,0.0950,0.1146,0.1319,0.1474,0.1532,0.1555];     %H04
    r_0 = [0.0160,0.0193,0.0230,0.0330,0.0440,0.0583,0.0728,0.0765,0.0785];     %N21

% when I used to do different shades:
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];   % chosen mix contour colors

map = [0.82 0.82 0.82; 1 1 1];
%map = [1 1 1];                     %when there isn't any solutions

map_neg = [0 0 0];

%range for post-Cr oxidation = modern day mantle FeO*
%0.02 to 0.06 contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.04 0.06 0.08 0.10 0.12 0.14 0.16];
%c0 = [0.02 0.03 0.04 0.05 0.06];
c0 = [0.02 0.06];
c_neg = [-0.05 0];                                  %for negative Fe3/sumFe limit
c_r0 = [-0.1, r_0(r0_val)-0.35/8.2];                %for original Fe3/sumFe value


figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, P(1:size(mix,1)), log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, P(1:size(nomix,1)), log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');

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
%text(5, -2.8, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
text(72, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
text(72, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)
% text(108, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
% text(118, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)

% annotation('textarrow',[0.6 0.57],[0.54 0.42])
% text(81, -1.6,'reduced MO','FontSize',8)
% annotation('textarrow',[0.61 0.64],[0.58 0.70])
% text(85, -1.3,'oxidized MO','FontSize',8)
title(ax1, title_name)


