% EL
% Aug 2022
% Updated 2024-02-22
%
% Efficiency vs Depth for last GI figure
%
% Note: need to change the data sheet because H04 and N21 have different
% ranges of pressure (line 18)
%
% also note, you may need to change the shaded region based on the max(mix)
% in line 68

reset = 0;      %reset=1 will read the file again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_Tconst.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'N21_data', 'Range', 1);
    P = data(2,:)/1e9;
    eff = data(1,1:28);
 
    H04_p5_nomix = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    H04_p5_mix = readmatrix(xlsx, 'Sheet', 'H04_5th_mix'); 
     
    H04_p25_nomix = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_p25_mix = readmatrix(xlsx, 'Sheet', 'H04_25th_mix'); 
    
    H04_p50_nomix = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_p50_mix = readmatrix(xlsx, 'Sheet', 'H04_50th_mix'); 
    
    H04_p75_nomix = readmatrix(xlsx, 'Sheet', 'H04_75th_nomix');
    H04_p75_mix = readmatrix(xlsx, 'Sheet', 'H04_75th_mix'); 

    N21_p5_nomix = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
    N21_p5_mix = readmatrix(xlsx, 'Sheet', 'N21_5th_mix'); 

    N21_p25_nomix = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_p25_mix = readmatrix(xlsx, 'Sheet', 'N21_25th_mix'); 

    N21_p50_nomix = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
    N21_p50_mix = readmatrix(xlsx, 'Sheet', 'N21_50th_mix'); 

    N21_p75_nomix = readmatrix(xlsx, 'Sheet', 'N21_75th_nomix');
    N21_p75_mix = readmatrix(xlsx, 'Sheet', 'N21_75th_mix');
    % 
end

mix = N21_p5_mix-0.35/8.2;       %0.35% reduction in FeO1.5 after Cr oxi with 8.2% FeO*
nomix = N21_p5_nomix-0.35/8.2;
r0_val = 3;
title_name = "N21: 5th Percentile";
letter = "d";

%Fe3/sumFe value BEFORE GI from modeling
%          [0th    1st    5th    25th   50th   75th   95th   99th   100th]
%Tconst method
    %r_0 = [0.0780,0.0828,0.0913,0.1043,0.1129,0.1211,0.1286,0.1311,0.1324];      %H04
    r_0 = [0.0382,0.0718,0.0834,0.0966,0.1059,0.1116,0.1239,0.1290,0.1310];      %N21
%Pmo method
    %r_0 = [0.0460,0.0549,0.0704,0.0987,0.1168,0.1327,0.1471,0.1515,0.1559];     %H04
    %r_0 = [0.0338,0.0383,0.0422,0.0602,0.0794,0.1050,0.1273,0.1346,0.1362];     %N21

% when I used to do different shades:
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];   % chosen mix contour colors

map = [0.82 0.82 0.82; 1 1 1];
%map = [1 1 1];                     %when there isn't any solutions

map_neg = [0 0 0];

%range for post-Cr oxidation = modern day mantle FeO*
%0.02 to 0.06 contours
c = linspace(0.01, 0.2, 20);
ct = [0.02 0.3 0.04 0.05 0.06 0.07 0.08 0.9 0.10 0.11 0.12];
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
%text(108, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')
%text(118, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)

% annotation('textarrow',[0.53 0.50],[0.34 0.23])
% text(59, -2.6,'reduced MO','FontSize',8)
% annotation('textarrow',[0.54 0.58],[0.36 0.48])
% text(69, -1.55,'oxidized MO','FontSize',8)
title(ax1, title_name)


