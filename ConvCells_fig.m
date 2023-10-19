% EL
% July 2022
% Updted 2023-10-06
% 
% Creates figure for convection cell mixing model, as fO2_surf vs time and
% Fe ratio vs time
% Uses data from CMBcells.xlsx, which is created by Ratios_CMBcells.m and
% data from MidMOcells.xlsx, which is created by Ratios_MidMOcells.m

reset = 1;

if reset == 1
    clear;

    % choose data rows
    row_10 = [2,3];
    row_20 = [4,5];
    row_40 = [6,7];
    row_P25 = [2,3];
    row_P50 = [4,5];
    row_P75 = [6,7];
    row_P100 = [8,9];

    data = readmatrix('CMBcells.xlsx', 'Sheet','3500');
    t = data(1,2:end);
    n10_r = data(row_10(1),2:end);
    n10_fO2 = data(row_10(2),2:end);
    n20_r = data(row_20(1),2:end);
    n20_fO2 = data(row_20(2),2:end);    
    n40_r = data(row_40(1),2:end);
    n40_fO2 = data(row_40(2),2:end);
    
    data = readmatrix('CMBcells.xlsx', 'Sheet','4500');
    n10_r_b = data(row_10(1),2:end);
    n10_fO2_b = data(row_10(2),2:end);
    n20_r_b = data(row_20(1),2:end);
    n20_fO2_b = data(row_20(2),2:end);
    n40_r_b = data(row_40(1),2:end);
    n40_fO2_b = data(row_40(2),2:end);
    
    data = readmatrix('MidMOcells.xlsx', 'Sheet','3500');
    P25_r = data(row_P25(1),2:end);
    P25_fO2 = data(row_P25(2),2:end);
    P50_r = data(row_P50(1),2:end);
    P50_fO2 = data(row_P50(2),2:end);
    P75_r = data(row_P75(1),2:end);
    P75_fO2 = data(row_P75(2),2:end);
    P100_r = data(row_P100(1),2:end);
    P100_fO2 = data(row_P100(2),2:end);
    
    data = readmatrix('MidMOcells.xlsx', 'Sheet','4500');
    P25_r_b = data(row_P25(1),2:end);
    P25_fO2_b = data(row_P25(2),2:end);
    P50_r_b = data(row_P50(1),2:end);
    P50_fO2_b = data(row_P50(2),2:end);
    P75_r_b = data(row_P75(1),2:end);
    P75_fO2_b = data(row_P75(2),2:end);
    P100_r_b = data(row_P100(1),2:end);
    P100_fO2_b = data(row_P100(2),2:end);

end

blue = [33,102,172;     %dark
    61, 121, 182;
    116, 159, 203;
    172, 198, 224;
    214, 227, 239]/255; %light
red = [213,62,79;
    218, 86, 101;
    229, 134, 145;
    239, 183, 189;
    248, 219, 222]/255;


figure(1);
subplot(1,2,1);     %Fe ratio
hold on
box on
plot(t/1e6, n10_r, 'Color', blue(1,:), 'LineStyle', '-', "LineWidth", 1.5);
plot(t/1e6, n20_r, 'Color', blue(1,:), 'LineStyle', '--', "LineWidth", 1.5);
plot(t/1e6, n40_r, 'Color', blue(1,:), 'LineStyle', ':', "LineWidth", 1.5);
plot(t/1e6, n10_r_b, 'Color', red(1,:), 'LineStyle', '-', "LineWidth", 1.5)
plot(t/1e6, n20_r_b, 'Color', red(1,:), 'LineStyle', '--', "LineWidth", 1.5)
plot(t/1e6, n40_r_b, 'Color', red(1,:), 'LineStyle', ':', "LineWidth", 1.5)
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.25])
yticks([0 0.05 0.10 0.15 0.20 0.25])
xlim([0 5])

subplot(1,2,2);     %fO2
hold on
box on
n10 = plot(t/1e6, n10_fO2, 'Color', blue(1,:), 'LineStyle', '-', "LineWidth", 1.5);
n20 = plot(t/1e6, n20_fO2, 'Color', blue(1,:), 'LineStyle', '--', "LineWidth", 1.5);
n40 = plot(t/1e6, n40_fO2, 'Color', blue(1,:), 'LineStyle', ':', "LineWidth", 1.5);
plot(t/1e6, n10_fO2_b, 'Color', red(1,:), 'LineStyle', '-', "LineWidth", 1.5)
plot(t/1e6, n20_fO2_b, 'Color', red(1,:), 'LineStyle', '--', "LineWidth", 1.5)
plot(t/1e6, n40_fO2_b, 'Color', red(1,:), 'LineStyle', ':', "LineWidth", 1.5)
xlabel("Time (Myr)")
ylabel('log(f_{O_2}) at surface')
ylim([-8 0])
yticks([-8 -6 -4 -2 0])
xlim([0 5])
%legend([n40 n20 n10], "n=10", "n=20", "n=40", "Location", "southeast")


figure(2);

subplot(2,2,1);     %3500K Fe ratio
hold on
box on
plot(t/1e6, P25_r, 'Color', blue(5,:), "LineWidth", 1.5);
plot(t/1e6, P50_r, 'Color', blue(4,:), "LineWidth", 1.5);
plot(t/1e6, P75_r, 'Color', blue(3,:), "LineWidth", 1.5);
plot(t/1e6, P100_r, 'Color', blue(2,:), "LineWidth", 1.5);
plot(t/1e6, n20_r, 'Color', blue(1,:), "LineWidth", 1.5);     %Pbase = 135
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.25])
yticks([0 0.05 0.10 0.15 0.20 0.25])
xlim([0 5])

subplot(2,2,2);     %3500K dIW
hold on 
box on
P25 = plot(t/1e6, P25_fO2, 'Color', blue(5,:), "LineWidth", 1.5);
P50 = plot(t/1e6, P50_fO2, 'Color', blue(4,:), "LineWidth", 1.5);
P80 = plot(t/1e6, P75_fO2, 'Color', blue(3,:), "LineWidth", 1.5);
P100 = plot(t/1e6, P100_fO2, 'Color', blue(2,:), "LineWidth", 1.5);
P135 = plot(t/1e6, n20_fO2, 'Color', blue(1,:), "LineWidth", 1.5);     %Pbase = 135
xlabel("Time (Myr)")
ylabel('log(f_{O_2}) at surface')
ylim([-8 0])
yticks([-8 -6 -4 -2 0])
xlim([0 5])
legend([P135 P100 P80 P50 P25], "135 GPa", "100 GPa", "80 GPa", "50 GPa", "25 GPa", "Location", "southeast")

subplot(2,2,3);     %4500K Fe ratio
hold on
box on
plot(t/1e6, P25_r_b, 'Color', red(5,:), "LineWidth", 1.5);
plot(t/1e6, P50_r_b, 'Color', red(4,:), "LineWidth", 1.5);
plot(t/1e6, P75_r_b, 'Color', red(3,:), "LineWidth", 1.5);
plot(t/1e6, P100_r_b, 'Color', red(2,:), "LineWidth", 1.5);
plot(t/1e6, n20_r_b, 'Color', red(1,:), "LineWidth", 1.5);     %Pbase = 135
xlabel("Time (Myr)")
ylabel('Fe^{3+}/\SigmaFe')
ylim([0 0.25])
yticks([0 0.05 0.10 0.15 0.20 0.25])
xlim([0 5])

subplot(2,2,4);     %4500K f
hold on
box on
P25 = plot(t/1e6, P25_fO2_b, 'Color', red(5,:), "LineWidth", 1.5);
P50 = plot(t/1e6, P50_fO2_b, 'Color', red(4,:), "LineWidth", 1.5);
P80 = plot(t/1e6, P75_fO2_b, 'Color', red(3,:), "LineWidth", 1.5);
P100 = plot(t/1e6, P100_fO2_b, 'Color', red(2,:), "LineWidth", 1.5);
P135 = plot(t/1e6, n20_fO2_b, 'Color', red(1,:), "LineWidth", 1.5);     %Pbase = 135
xlabel("Time (Myr)")
ylabel('log(f_{O_2}) at surface')
ylim([-8 0])
yticks([-8 -6 -4 -2 0])
xlim([0 5])
legend([P135 P100 P80 P50 P25], "135 GPa", "100 GPa", "80 GPa", "50 GPa", "25 GPa", "Location", "southeast")

