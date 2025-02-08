% EL
% Aug 2022
% Updated 2024-07-03
%
% Efficiency vs Depth for last GI figure. Data from EffvsD.m. Creates
% Figure 3 in MO redox manuscript (and the Supp Fig Pmo version)
%
% Change the following:
%   xlsx        data file, either **_const or **_Pmo
%   data        data sheet because H04 and N21 have different ranges of pressure
%   mix & nomix which data to use
%   r0_val      which initial Fe3+/sumFe value from pre-GI MC results,
%                   comment out the r0 arrays you aren't using.
%   title       for graph
%   letter      for graph
%
% also note, you may need to change the shaded contour region based on the max(mix)
% and the location of the letter and 'r0=##'. And to be honest, I often
% just manually move around the oxidized/reduced arrows.


reset = 0;      %reset=1 will read the file again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_const.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'N21_data', 'Range', 1);
    P = data(2,:)/1e9;
    eff = data(1,1:28);

    H04_p1_nomix = readmatrix(xlsx, 'Sheet', 'H04_1st_nomix');
    H04_p1_mix = readmatrix(xlsx, 'Sheet', 'H04_1st_mix'); 
 
    H04_p5_nomix = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    H04_p5_mix = readmatrix(xlsx, 'Sheet', 'H04_5th_mix'); 
     
    H04_p25_nomix = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_p25_mix = readmatrix(xlsx, 'Sheet', 'H04_25th_mix'); 
    
    H04_p50_nomix = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_p50_mix = readmatrix(xlsx, 'Sheet', 'H04_50th_mix'); 
    
    % H04_p75_nomix = readmatrix(xlsx, 'Sheet', 'H04_75th_nomix');
    % H04_p75_mix = readmatrix(xlsx, 'Sheet', 'H04_75th_mix');

    N21_p1_nomix = readmatrix(xlsx, 'Sheet', 'N21_1st_nomix');
    N21_p1_mix = readmatrix(xlsx, 'Sheet', 'N21_1st_mix'); 

    N21_p5_nomix = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
    N21_p5_mix = readmatrix(xlsx, 'Sheet', 'N21_5th_mix'); 

    N21_p25_nomix = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_p25_mix = readmatrix(xlsx, 'Sheet', 'N21_25th_mix'); 

    N21_p50_nomix = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
    N21_p50_mix = readmatrix(xlsx, 'Sheet', 'N21_50th_mix'); 

    % N21_p75_nomix = readmatrix(xlsx, 'Sheet', 'N21_75th_nomix');
    % N21_p75_mix = readmatrix(xlsx, 'Sheet', 'N21_75th_mix');
    % 
end

mix = N21_p5_mix-0.35/8.2;       %0.35% reduction in FeO1.5 after Cr oxi with 8.2% FeO*
nomix = N21_p5_nomix-0.35/8.2;
r0_val = 3;
title_name = "N21: 5th Percentile";
letter = "d";

%CHOOSE YOUR Fe3/sumFe VALUE BEFORE GI
%         [0th    1st    5th    25th   50th   75th   95th   99th   100th]
%H04
    %r_0 = [0.0695,0.0781,0.0878,0.1009,0.1101,0.1189,0.1271,0.1299,0.1314];    %Tconst
    %r_0 = [0.0365,0.0523,0.0709,0.0963,0.1156,0.1316,0.1470,0.1510,0.1538];    %Pmo
%N21
    r_0 = [0.0569,0.0644,0.0804,0.0953,0.1047,0.1108,0.1210,0.1267,0.1289];   %Tconst
    %r_0 = [0.0311,0.0367,0.0426,0.0615,0.0793,0.1052,0.1285,0.1338,0.1362];    %Pmo

% when I used to do different shades:
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];   % chosen mix contour colors

map = [0.82 0.82 0.82; 1 1 1];
%map = [1 1 1];                     %when there isn't any solutions

map_neg = [0 0 0];

%range for post-Cr oxidation = modern day mantle FeO*
%0.02 to 0.06 contours
c = linspace(0.01, 0.2, 20);
%ct = [0.02 0.3 0.04 0.05 0.06 0.07 0.08 0.9 0.10 0.11 0.12];
ct = [0.02 0.04 0.06 0.08 0.10 0.12];               %only even contours labeled
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
text(96, -0.15, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')    %N21
text(96, -0.35, letter, 'FontWeight', 'bold', 'FontSize', 20)           %N21
%text(112, -2.88, "r_0=" + round(r_0(r0_val),3), 'FontWeight', 'bold')   %H04
%text(122, -2.65, letter, 'FontWeight', 'bold', 'FontSize', 20)          %H04

% annotation('textarrow',[0.56 0.52],[0.45 0.30])
% text(64.2, -2.38,'reduced MO','FontSize',8)
% annotation('textarrow',[0.56 0.60],[0.47 0.61])
% text(71.5, -1.1,'oxidized MO','FontSize',8)
title(ax1, title_name)


