% ELH
% 02-15-2024

% parameterizing Rubie+2011 Fig 4b and using data from Table S3a to get the
% impactor Fe3+/sumFe and Earth composition over time (tracking FeO
% content).

% outputs a file for impactor SiO2, FeO, Fe3+/sumFe, and Earth composition
% over time for a given model.

clear;

% INPUTS
model = 4;                      %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
write = 1;
    fileOut = "Compositions.xlsx";
    sheetNameE = "H04_E";
    sheetNameImp = "H04_Imp";

% Rubie+2011 FeO data for accreting Earth, from Table S3a
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11');
F_temp = data(13, 4:end);      % earth mass fraction
OxiWts_temp = data(1:12,4:end);
OxiMolW = data(1:12,3);     %  note need the full molecular weight here, not the 'per oxide'

% dIW_R11 = [-5.09, -5.25, -5.27, -5.27, -5.18, -5.09, -4.96, -4.85, -4.73, -4.59, ...
%     -4.32, -4.18, -4.00, -3.85, -3.70, -3.55, -3.40, -3.24, -3.12, -3.00, -2.60, ...
%     -2.35, -2.15, -1.98];

% SET UP ACCRETION MODEL
M_E_0 = 5.97e24;                                %[kg] present day mass of Earth
M_c_0 = 1.88e24;                                %[kg] present day mass of core
[t,Accr_model] = getAccrModel(model);           %Accr_model = F = mass fraction
M_E = M_E_0 * Accr_model;
M_imp = M_E_0 * diff(Accr_model);

OxiWts = zeros(length(OxiMolW), length(Accr_model));
for i = 1:length(Accr_model)
    if Accr_model(i) < 0.11
        OxiWts(5,i) = 0.3/0.11 * Accr_model(i);       %linear initially to 'extrap' from Rubie+11 when F<0.11
        OxiWts(1,i) = 100 - OxiWts(5,i) - sum(OxiWts_temp(2:4,1)) - sum(OxiWts_temp(6:12,1)); %compensate FeO drop with Si increase
        OxiWts(2:4,i) = OxiWts_temp(2:4,1);
        OxiWts(6:12,i) = OxiWts_temp(6:12,1);
    elseif Accr_model(i) < 0.67
        for j = 1:length(OxiMolW)
            OxiWts(j,i) = interp1(F_temp(1:20), OxiWts_temp(j,1:20), Accr_model(i), 'pchip');
        end
    else
        for j = 1:length(OxiMolW)
            OxiWts(j,i) = interp1(F_temp(20:24), OxiWts_temp(j,20:24), Accr_model(i), 'pchip');
        end
    end
end


% figure(1);      %check if Accr_model values follow Rubie+11 FeO
% hold on
% box on
% plot(Accr_model, OxiWts, 'MarkerSize', 8)
% xlim([0 1])
% xlabel("mass fraction")
% ylabel("wt %")

% check that this reproduces dIW okay with changing FeO content like this.
% Find mole fraction of FeO for Earth based on varying composition
tempMol = OxiWts./OxiMolW;
MolSum = sum(tempMol,1);    %total moles
Mol_E = tempMol./MolSum;     %mole fraction
FeO_Mol_E = Mol_E(5,:);      %extract FeO mole fraction
Fe_Mol = 0.78;              %estimate mole fraction for Fe in metal

dIW_E_simple = 2*log10(FeO_Mol_E/Fe_Mol);
dIW_E_simple(isinf(dIW_E_simple)) = 0;
dIW_E_long = 2*log10((1.148*FeO_Mol_E + 1.319*FeO_Mol_E.^2)/Fe_Mol);     %with Rubie+11 param to X for wustite
dIW_E_long(isinf(dIW_E_long)) = 0;

figure(2);
hold on
box on
plot(Accr_model, dIW_E_simple, 'b.', "MarkerSize", 8, "DisplayName", "simple model")
plot(Accr_model, dIW_E_long, 'r.', "MarkerSize", 8, "DisplayName", "with R11 params")
plot(F_temp, dIW_R11, 'k.', "MarkerSize", 8, "DisplayName", "from R11 graph")
xlabel("Mass fraction")
ylabel("dIW")
title("Earth")
legend("Location", "southeast")

% assume core and mantle take up proportional mass of Earth
M_c = M_E * M_c_0/M_E_0;    %core
M_m = M_E - M_c;            %mantle
Fe_ratio = 0.321;               % [] weight % of iron on Earth/impactor
Si_ratio = 1-Fe_ratio;          % [] weight % of silicate 'mantle' on Earth/impactor
M_m_imp = M_imp*Si_ratio;       % [kg] approximate proportion silicate mass of impactor

FeO_tot = M_m.*OxiWts(5,:)/100;           %[kg] mass of FeO for Earth
FeO_added = diff(FeO_tot);      %[kg] difference of FeO mass between impacts

FeO_imp = FeO_added./M_m_imp*100;             %[ratio] FeO content of the impactor mantle
FeO_imp(isnan(FeO_imp)) = 0;

% Find mole fraction of FeO for impactor based on varying composition
% (changes in SiO2 compensated by changes in FeO)
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpLate');
OxiWts_late_temp = data(:,2);
OxiMolW = data(:,4);
OxiMolW_byM = data(:,3);
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpEarly');
OxiWts_early_temp = data(:,2);

Si_imp = zeros(1, length(FeO_imp));
FeO_Mol_imp = zeros(1, length(FeO_imp));
OxiWts_late = zeros(length(OxiWts_late_temp), length(FeO_imp))+OxiWts_late_temp;    %keep other comps the same, but change Si and Fe
OxiWts_early = zeros(length(OxiWts_early_temp), length(FeO_imp))+OxiWts_early_temp;
for i = 1:length(FeO_imp)
    if FeO_imp(i) > 0 && Accr_model(i+1) < 0.67
        OxiWts_early(1,i) = 100 - (OxiWts_early(3,i)+OxiWts_early(7,i)+OxiWts_early(9,i)+FeO_imp(i));  %SiO2 content
        Si_imp(i) = OxiWts_early(1,i);
        OxiWts_early(5,i) = FeO_imp(i);                                                                %FeO wt% content
        tempMol = OxiWts_early(:,i)./OxiMolW;
        MolSum = sum(tempMol);    %total moles
        Mol = tempMol/MolSum;     %mole fraction
        FeO_Mol_imp(i) = Mol(5);      %extract FeO mole fraction
    elseif FeO_imp(i) > 0 && Accr_model(i+1) >= 0.67
        OxiWts_late(1,i) = 100 - (OxiWts_late(3,i)+OxiWts_late(7,i)+OxiWts_late(9,i)+FeO_imp(i));  %SiO2 content
        Si_imp(i) = OxiWts_late(1,i);
        OxiWts_late(5,i) = FeO_imp(i);                                                             %FeO content
        tempMol = OxiWts_late(:,i)./OxiMolW;
        MolSum = sum(tempMol);
        Mol = tempMol/MolSum;
        FeO_Mol_imp(i) = Mol(5);
    else
        FeO_Mol_imp(i) = 0;   %so that we only have impactor times
    end
end

dIW_imp_simple = 2*log10(FeO_Mol_imp/Fe_Mol);
dIW_imp_simple(isinf(dIW_imp_simple)) = 0;
dIW_imp_long = 2*log10((1.148*FeO_Mol_imp + 1.319*FeO_Mol_imp.^2)/Fe_Mol);     %with Rubie+11 param to X for wustite
dIW_imp_long(isinf(dIW_imp_long)) = 0;

% figure(3);
% hold on
% box on
% plot(Accr_model(2:end), dIW_imp_simple, 'b.', "MarkerSize", 8, "DisplayName", "simple model")
% plot(Accr_model(2:end), dIW_imp_long, 'r.', "MarkerSize", 8, "DisplayName", "with R11 params")
% xlabel("Mass fraction")
% ylabel("dIW")
% title("Impactor")
% legend("Location", "southeast")
% legend

%mean values from impactor T and P equilibrium of Table S3a, Rubie+2011
T_early = 2226;
P_early = 5.4;
T_late = 2324;
P_late = 8.8;

PV_table = readmatrix('\db\PVcalc_old.xlsx');
%PV_term = calcPV(T_early, P_early, PV_table);      %need to see if I can use this instead
P_test = PV_table(:,1);
T_test = [2000, 2500, 3000, 3500, 4000, 4500];
PV_test = PV_table(:,2:end);

PV_temp = interp2(T_test,P_test,PV_test,T_early,P_early);
PV_term_early = PV_temp/log(10);
PV_temp = interp2(T_test,P_test,PV_test,T_late,P_late);
PV_term_late = PV_temp/log(10);

FeRatio_imp = zeros(1,length(FeO_imp));              % determine Fe3+/sumFe of impactor based on P, T, dIW
for i = 1:length(FeO_imp)
    if FeO_imp(i) > 0 && Accr_model(i+1) < 0.67
        FeRatio_imp(i) = calcFeRatio_byIW(T_early, P_early, dIW_imp_long(i), PV_term_early, OxiWts_early(:,i), OxiMolW_byM);
    elseif FeO_imp(i) > 0 && Accr_model(i+1) >= 0.67
        FeRatio_imp(i) = calcFeRatio_byIW(T_late, P_late, dIW_imp_long(i), PV_term_late, OxiWts_late(:,i), OxiMolW_byM);
    else
        FeRatio_imp(i) = 0;
    end
end


if write == 1
    writematrix(Si_imp, fileOut, 'Sheet', sheetNameImp)
    writematrix(FeO_imp, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')
    writematrix(FeRatio_imp, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')
    writematrix(OxiWts, fileOut, 'Sheet', sheetNameE)
end

FeRatio_imp(FeRatio_imp==0) = NaN;   % so it doesn't plot the zeros

figure(4);
hold on
box on
plot(t(2:end), FeRatio_imp, 'b.', "MarkerSize", 8)
ylabel("Fe3+/sumFe")
xlabel("Time (Myr)")
