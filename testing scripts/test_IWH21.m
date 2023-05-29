% EL
% Feb 2023
%
% Testinng H21 IW buffer past 130 GPa
% Compare equation evaluation and point linear extrapolation
% Hirschmann states that his parameterization is "not recommended" beyond
% 100 GPa

clear;

data = readmatrix('/db/geotherms.xlsx', 'Sheet', '3000K');
T = data(:,4);          % to 135 GPa
T_short = T(1:226);     % truncated @ 100 GPa
P = data(:,2);
P_short = P(1:226);

T_const = 3000;

% Evaluated IW using Hirschmann 2021 IW function
IW_eval = zeros(length(P),1);
for i = 1:length(P)
    IW_eval(i) = getIW_H21(P(i),T(i));
end

IW_eval_c = zeros(length(P),1);
for i = 1:length(P)
    IW_eval_c(i) = getIW_H21(P(i),T_const);
end

% Linearly extrapolated >100 GPa
IW_lin = interp1(P_short, IW_eval(1:226), P(227:end), 'linear', 'extrap');
IW_lin_c = interp1(P_short, IW_eval_c(1:226), P(227:end), 'linear', 'extrap');

%Campbell et al. 2009;
a0 = 6.54106; a1 = 0.0012324;
b0 = -28163.6; b1 = 546.32; b2 = -1.13412; b3 = 0.0019274;
IW_C09 = (a0+a1*P) + (b0+b1*P+b2*(P.^2)+b3*(P.^3))./T;

a0 = 6.54106; a1 = 0.0012324;
b0 = -28163.6; b1 = 546.32; b2 = -1.13412; b3 = 0.0019274;
IW_C09_c = (a0+a1*P) + (b0+b1*P+b2*(P.^2)+b3*(P.^3))./T_const;


IW_diff = IW_lin - IW_eval(227:end);

figure(1);
hold on
box on
plot(P, IW_eval, 'LineWidth', 1.5, 'DisplayName', 'H21 Evaluated (Geo)')
plot(P(227:end), IW_lin, 'LineWidth', 1.5, 'DisplayName', 'H21 Extrapolated (Geo)')
plot(P, IW_C09, 'LineWidth', 1.5, 'DisplayName', 'Campbell+09 (Geo)')
plot(P, IW_eval_c, ':', 'Color', '#0072BD', 'LineWidth', 1.5, 'DisplayName', 'H21 Evaluated (Const)')
plot(P(227:end), IW_lin_c, ':', 'Color', '#D95319', 'LineWidth', 1.5, 'DisplayName', 'H21 Extrapolated (Const)')
plot(P, IW_C09_c, ':', 'Color', '#EDB120', 'LineWidth', 1.5, 'DisplayName', 'Campbell+09 (Const)')
xlim([50, 140])
xlabel('Pressure (GPa)')
ylabel('log(IW)')
legend('Location', 'southeast')
hold off