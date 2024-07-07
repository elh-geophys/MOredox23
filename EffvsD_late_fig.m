% EL
% Aug 2022
% Updated 2023-07-03
%
% Efficiency vs Depth figure (Late Accretion)
%
% Adjust the following:
%   xlsx            which spredsheet to use
%   mix & nomix     which data sheet to use
%   r_0_idx         which Fe3+/sumFe ratio was the initial one, comment out
%                       the array not using
%   title           for graph
%   letter          for graph
%
% Note you may need to adjust the contour shading if there are no solutions
% (check min(mix)... if nothing less than 0.6, then be sad :( ). Lastly, I
% usually manually move the reduced/oxidized arrows instead of worrying to
% much about scripting their exact location.

reset = 0;     %if data sheets need to be read again

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_late_Tconst.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'data', 'Range', 1);
    P = data(2,1:100)/1e9;      %just to ~50GPa
    eff = data(1,10:28);        %just 1%-100%            

    %these are all set before Cr oxi    
    data = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_50th_nomix = data(1:100,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_50th_mix');
    H04_50th_mix = data(1:100,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_25th_nomix = data(1:100,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_25th_mix');
    H04_25th_mix = data(1:100,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    H04_5th_nomix = data(1:100,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_5th_mix');
    H04_5th_mix = data(1:100,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'H04_1st_nomix');
    H04_1st_nomix = data(1:100,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_1st_mix');
    H04_1st_mix = data(1:100,10:28);
    
%     data = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
%     N21_50th_nomix = data(1:100,10:28);
%     data = readmatrix(xlsx, 'Sheet', 'N21_50th_mix');
%     N21_50th_mix = data(1:100,10:28);
%     
%     data = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
%     N21_25th_nomix = data(1:100,10:28);
%     data = readmatrix(xlsx, 'Sheet', 'N21_25th_mix');
%     N21_25th_mix = data(1:100,10:28);
%     
%     data = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
%     N21_5th_nomix = data(1:100,10:28);
%     data = readmatrix(xlsx, 'Sheet', 'N21_5th_mix');
%     N21_5th_mix = data(1:100,10:28);
    
end

mix = H04_1st_mix - 0.35/8.2;       %0.35% reduction in FeO1.5 after Cr oxi with 8.1% FeO*
nomix = H04_1st_nomix - 0.35/8.2;
r_0_idx = 1;                            %for r_0 to choose
title_name = "H04: 1st Percentile";
letter = "a";

%Fe3/sumFe value AFTER GI
%Tconst method
         % 1st     5th     25th    50th
r_0_temp = [0.0856,0.0955,0.1107,0.1200];    %Tconst
%r_0_temp = [0.0613,0.0824,0.1077,0.1260];    %Pmo

r_0 = r_0_temp(r_0_idx)-0.35/8.2;

% back when I used to do shades of gray:
%       0.2-0.3         0.3-0.4         0.4-0.5         0.5-0.6         rest
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
        
map = [0.82 0.82 0.82; 1 1 1];          % chosen mix contour colors
%map = [1 1 1; 0.82 0.82 0.82];
%map = [1 1 1];                           % for no solution, check min(mix)

map_neg = [0 0 0];

c = linspace(0.005, 0.1, 20);
ct = [0.01 0.02 0.03 0.04 0.05 0.06];
%c0 = [0.02 0.03 0.04 0.05 0.06];
c0 = [0.02 0.06];
c_neg = [-0.05 0];                      %for negative Fe3/sumFe limit
c_r0 = [-0.1, r_0];                     %for original Fe3/sumFe value

figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, P, log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, P, log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 500, 'FaceColor', 'none');

ax3 = axes;
contour(ax3, P, log10(eff), nomix', c_neg, 'LineWidth', 2, 'LineStyle', ':', 'FaceColor', 'none');

ax4 = axes;
contour(ax4, P, log10(eff), nomix', c_r0, 'LineWidth', 1, 'LineStyle', '--', 'Color', 'b', 'FaceColor', 'none');

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
ax1.XLim = [min(P) 50];
ax1.XLabel.String = 'Pressure (GPa)';
ax1.YLabel.String = 'Degree of Chemical Equilibration';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};
text(42, -1.87, "r_0=" + round(r_0,3), 'FontWeight', 'bold')
text(47, -1.70, letter, 'FontWeight', 'bold', 'FontSize', 20)
% annotation('textarrow',[0.67 0.59],[0.51 0.39])
% text(24, -1.37,'reduced MO','FontSize',8)
% annotation('textarrow',[0.68 0.75],[0.53 0.65])
% text(34, -0.6,'oxidized MO','FontSize',8)
title(ax1, title_name)

