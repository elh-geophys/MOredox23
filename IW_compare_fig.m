% EL
% Feb 2023
%
% Testinng H21 IW buffer past 130 GPa
% Compare equation evaluation and point linear extrapolation
% Hirschmann states that his parameterization is "not recommended" beyond
% 100 GPa

clear;

geotherms = ["2500K", "3500K", "4500K"];

T = zeros(290,length(geotherms));          % to 135 GPa
P = zeros(290,length(geotherms));
T_short = zeros(226,length(geotherms));    % truncated at 100 GPa
P_short = zeros(226,length(geotherms));

for j = 1:length(geotherms)
    data = readmatrix('/db/geotherms.xlsx', 'Sheet', geotherms(j));
    T(:,j) = data(2:end,4);         
    T_short(:,j) = T(1:226);     
    P(:,j) = data(2:end,2);
    P_short(:,j) = P(1:226);
end

% Evaluated IW using Hirschmann 2021 IW function
IW_eval = zeros(length(P),length(geotherms));
for j = 1:length(geotherms)                               %j index for geotherms, i index for P
    for i = 1:length(P)
        IW_eval(i,j) = getIW_H21(P(i,j)*1e9,T(i,j));
    end
end

% Linearly extrapolated >100 GPa
IW_lin = zeros(length(P(227:end,1)),length(geotherms));
for j = 1:length(geotherms)
    IW_lin(:,j) = interp1(P_short(:,j), IW_eval(1:226,j), P(227:end,j), 'linear', 'extrap');
end

% Campbell et al. 2009
a0 = 6.54106; a1 = 0.0012324;
b0 = -28163.6; b1 = 546.32; b2 = -1.13412; b3 = 0.0019274;
IW_C09 = zeros(length(P), length(geotherms));
for j = 1:length(geotherms)
    IW_C09(:,j) = (a0+a1*P(:,j)) + (b0+b1*P(:,j)+b2*(P(:,j).^2)+b3*(P(:,j).^3))./T(:,j);
end

%IW_diff = IW_lin - IW_eval(227:end);

figure(1);
hold on
box on
colors = ["#0072BD", "#EDB120", "#77AC30"];
ref = zeros(length(geotherms),1);
for j = 1:length(geotherms)
    ref(j) = plot(P(:,j), IW_eval(:,j), '-', 'Color', colors(j), 'LineWidth', 1.5);
    plot(P(227:end,j), IW_lin(:,j), ':', 'Color', colors(j), 'LineWidth', 1.5)
    plot(P(:,j), IW_C09(:,j), '--', 'Color', colors(j), 'LineWidth', 1.5)
end
xlim([90, 135])
xticks([90, 100, 110, 120, 130])
yticks([9, 10, 11, 12, 13])
xlabel('Pressure (GPa)')
ylabel('log(IW)')
legend(ref, "2500K", "3500K", "4500K", 'Location', 'southeast')
hold off