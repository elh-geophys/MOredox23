% EL
% Aug 2022
% Updated 2023-10-02
%
% Efficiency vs Depth figure (Late Accretion)

reset = 0;

if reset == 1
    clear;
    
    xlsx = '\db\Rain_EffvsD_late_U2Q.xlsx';
    
    data = readmatrix(xlsx, 'Sheet', 'data', 'Range', 1);
    P = data(2,1:97)/1e9;      %just to ~50GPa
    eff = data(1,10:28);        %just 1%-100%            

    %these are all set before Cr oxi    
    data = readmatrix(xlsx, 'Sheet', 'H04_50th_nomix');
    H04_50th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_50th_mix');
    H04_50th_mix = data(:,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'H04_25th_nomix');
    H04_25th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_25th_mix');
    H04_25th_mix = data(:,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'H04_5th_nomix');
    H04_5th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'H04_5th_mix');
    H04_5th_mix = data(:,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'N21_50th_nomix');
    N21_50th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'N21_50th_mix');
    N21_50th_mix = data(:,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'N21_25th_nomix');
    N21_25th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'N21_25th_mix');
    N21_25th_mix = data(:,10:28);
    
    data = readmatrix(xlsx, 'Sheet', 'N21_5th_nomix');
    N21_5th_nomix = data(:,10:28);
    data = readmatrix(xlsx, 'Sheet', 'N21_5th_mix');
    N21_5th_mix = data(:,10:28);
    
end

mix = H04_5th_mix - 0.35/8.1;       %0.35% reduction in FeO1.5 after Cr oxi with 8.1% FeO* from Deng20 comp
nomix = H04_5th_nomix - 0.35/8.1;
r_0_idx = 1;                            %for r_0 to choose
title_name = "H04: 5th Percentile";
letter = "a";

%Fe3/sumFe value AFTER GI
          % 5th     25th    50th
%U2Q method
    r_0 = [0.0507, 0.0641, 0.0808];      %H04
    %r_0 = [0.0901, 0.1033, 0.1189];      %N21
    
r_0 = r_0(r_0_idx)-0.35/8.1;
     
%map = [0.88 0.88 0.88; 0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
%map = [0.88 0.88 0.88; 0.82 0.82 0.82];
%map = [0.88 0.88 0.88; 0.88 0.88 0.88; 1 1 1];
%map = [0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72; 1 1 1];
map = [0.72 0.72 0.72; 1 1 1; 1 1 1];
%map = [0.82 0.82 0.82; 0.78 0.78 0.78; 0.72 0.72 0.72];  
%map = [0.78 0.78 0.78; 0.72 0.72 0.72];                % chosen mix contour colors

map_neg = [0 0 0];

c = linspace(0.005, 0.1, 20);
ct = [0.01 0.02 0.03 0.04 0.05 0.06];
%c0 = [0.02 0.03 0.04 0.05 0.06];
c0 = [0.02 0.03 0.04 0.05 0.06];
c_neg = [-0.03 0];

figure('Position', [200 200 500 400]);

ax1 = axes;
contourf(ax1, P(1:size(mix,1)), log10(eff), mix', c0, 'LineWidth', 2, 'EdgeColor', 'none');

ax2 = axes;
contour(ax2, P(1:size(mix,1)), log10(eff), nomix', c, 'LineWidth', 2, 'ShowText', 'on', 'TextList', ct, 'LabelSpacing', 200, 'FaceColor', 'none');

ax3 = axes;
contour(ax3, P(1:size(mix,1)), log10(eff), nomix', c_neg, 'LineWidth', 1, 'LineStyle', '--', 'FaceColor', 'none');

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
ax1.XLim = [min(P) P(size(mix,1))];
ax1.XLabel.String = 'Pressure (GPa)';
ax1.YLabel.String = 'Degree of Chemical Equilibration';
ax1.YTick = log10(eff);
ax1.YTickLabel = {'1%' '' '' '' '' '' '' '' '' ...
    '10%' '' '' '' '' '' '' '' '' ...
    '100%'};
text(14.4, -1.87, "r_0=" + round(r_0,3), 'FontWeight', 'bold')
text(15.9, -1.70, letter, 'FontWeight', 'bold', 'FontSize', 20)
title(ax1, title_name)

