% ELH
% 02-15-2024
% Updated 6/30/24

% parameterizing Rubie+2011 Fig 4b and using data from Table S3a and S3b
% to get the impactor Fe3+/sumFe and Earth composition over time (tracking FeO
% content).

% outputs a file for impactor SiO2, FeO, Fe3+/sumFe, and Earth composition
% over time for a given model.

clear;
clf;

% INPUTS
model = 4;                      %accretion models, 4 = H04, 5 = N21
write = 0;
    fileOut = "Compositions.xlsx";
    sheetNameE = "H04_E";
    sheetNameImp = "H04_Imp";

% Rubie+2011 FeO data for accreting Earth, from Table S3a
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Emantle');
F_temp = data(13, 4:end);      % earth mass fraction
OxiWts_temp = data(1:12,4:end);
OxiMolWts_byM = data(1:12,2);
OxiMolWts = data(1:12,3);     %  note need the full molecular weight here, not the 'per oxide'

% Rubie+2011 core data for accreting Earth, from Table S3b
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'Rubie11_Ecore');
MetalMolWts = data(1:6,2);
MetalWts_temp = data(1:6, 3:end);
M_c_ratio_temp = data(8,3:end);

dIW_R11 = [-5.09, -5.25, -5.27, -5.27, -5.18, -5.09, -4.96, -4.85, -4.73, -4.59, ...
     -4.32, -4.18, -4.00, -3.85, -3.70, -3.55, -3.40, -3.24, -3.12, -3.00, -2.60, ...
     -2.35, -2.15, -1.98];

% SET UP ACCRETION MODEL
M_E_0 = 5.97e24;                                %[kg] present day mass of Earth
[t,Accr_model] = getAccrModel(model);           %Accr_model = F = mass fraction
M_E = M_E_0 * Accr_model;
M_imp = M_E_0 * diff(Accr_model);

OxiWts = zeros(length(OxiMolWts), length(Accr_model));
MetalWts = zeros(length(MetalMolWts), length(Accr_model));
E2L_idx = 0;    % to track the idx for the early-to-late transition for composition
preGI_idx = 0;  % to track the idx for pre-GI (for N21 model)

if model == 5
    Accr_model = 2*Accr_model;    %for N21, last GI is "Earth-like" so scheme is kinda 'scrunched' at the initial part of accretion
    for i = 1:length(Accr_model)
        if Accr_model(i) < 0.11   %the N21 doesn't actually need this, just here for completionist sake
            OxiWts(5,i) = 0.3/0.11 * Accr_model(i);       %linear initially to 'extrap' from Rubie+11 when F<0.11
            OxiWts(1,i) = 100 - OxiWts(5,i) - sum(OxiWts_temp(2:4,1)) - sum(OxiWts_temp(6:12,1)); %compensate FeO drop with Si increase
            OxiWts(2:4,i) = OxiWts_temp(2:4,1);
            OxiWts(6:12,i) = OxiWts_temp(6:12,1);

            MetalWts(:,i) = MetalWts_temp(:,1);    %keep core comp the same as initial since ~constant early on
            
            E2L_idx = E2L_idx + 1;
        elseif Accr_model(i) < 0.67     %early impactors
            for j = 1:length(OxiMolWts)
                OxiWts(j,i) = interp1(F_temp(1:20), OxiWts_temp(j,1:20), Accr_model(i), 'pchip');
            end

            for k = 1:length(MetalMolWts)
                MetalWts(k,i) = interp1(F_temp(1:20), MetalWts_temp(k,1:20), Accr_model(i), 'pchip');
            end

            E2L_idx = E2L_idx + 1;
        elseif Accr_model(i) <= 1.0     %late impactors
            for j = 1:length(OxiMolWts)
                OxiWts(j,i) = interp1(F_temp(20:24), OxiWts_temp(j,20:24), Accr_model(i), 'pchip');
            end

            for k = 1:length(MetalMolWts)
                MetalWts(k,i) = interp1(F_temp(20:24), MetalWts_temp(k,20:24), Accr_model(i), 'pchip');
            end
            preGI_idx = preGI_idx + 1;
        else                           %GI "Earth-like"                                                       
            for j = 1:length(OxiMolWts)
                OxiWts(j,i) = interp1(F_temp(20:24), OxiWts_temp(j,20:24), 0.99, 'pchip');
            end

            for k = 1:length(MetalMolWts)
                MetalWts(k,i) = interp1(F_temp(20:24), MetalWts_temp(k,20:24), 0.99, 'pchip');
            end
        end
    end
    Accr_model = 0.5*Accr_model;
    preGI_idx = preGI_idx + E2L_idx;
else
    for i = 1:length(Accr_model)
        if Accr_model(i) < 0.11
            OxiWts(5,i) = 0.3/0.11 * Accr_model(i);       %linear initially to 'extrap' from Rubie+11 when F<0.11
            OxiWts(1,i) = 100 - OxiWts(5,i) - sum(OxiWts_temp(2:4,1)) - sum(OxiWts_temp(6:12,1)); %compensate FeO drop with Si increase
            OxiWts(2:4,i) = OxiWts_temp(2:4,1);
            OxiWts(6:12,i) = OxiWts_temp(6:12,1);

            MetalWts(:,i) = MetalWts_temp(:,1);    %keep core comp the same as initial since ~constant early on
            
            E2L_idx = E2L_idx + 1;
        elseif Accr_model(i) < 0.67     %early impactors
            for j = 1:length(OxiMolWts)
                OxiWts(j,i) = interp1(F_temp(1:20), OxiWts_temp(j,1:20), Accr_model(i), 'pchip');
            end

            for k = 1:length(MetalMolWts)
                MetalWts(k,i) = interp1(F_temp(1:20), MetalWts_temp(k,1:20), Accr_model(i), 'pchip');
            end

            E2L_idx = E2L_idx + 1;
        else                            %late impactors, incl GI
            for j = 1:length(OxiMolWts)
                OxiWts(j,i) = interp1(F_temp(20:24), OxiWts_temp(j,20:24), Accr_model(i), 'pchip');
            end

            for k = 1:length(MetalMolWts)
                MetalWts(k,i) = interp1(F_temp(20:24), MetalWts_temp(k,20:24), Accr_model(i), 'pchip');
            end
        end
    end
end

% NOTE: make F_temp/2 for N21 checks!!!
figure(1);      %check if Accr_model values follow Rubie+11 FeO and Fe
subplot(1,2,1)
hold on
box on
plot(Accr_model, OxiWts(5,:), 'r.', 'MarkerSize', 8, 'DisplayName', 'Calculated')
plot(F_temp, OxiWts_temp(5,:), 'k-', 'DisplayName', 'From R11 data')
xlim([0 1])
xlabel("Mass Fraction")
ylabel("Wt %")
title("FeO in Earth Mantle vs Accretion")

subplot(1,2,2)
hold on
box on
plot(Accr_model, MetalWts(1,:), 'r.', 'MarkerSize', 8)
plot(F_temp, MetalWts_temp(1,:), 'k-')
xlim([0 1])
xlabel("Mass Fraction")
ylabel("Wt %")
title("Fe in Earth Core vs Accretion")

% check that this reproduces dIW okay with changing FeO content like this.
% Find mole fraction of FeO for Earth based on varying composition
tempMol = OxiWts./OxiMolWts;
MolSum = sum(tempMol,1);    %total moles
Mol_E = tempMol./MolSum;     %mole fraction
FeO_Mol_E = Mol_E(5,:);      %extract FeO mole fraction

tempMol = MetalWts./MetalMolWts;
MolSum = sum(tempMol,1);
Mol_Ecore = tempMol./MolSum;
Fe_Mol_core = Mol_Ecore(1,:);              %estimate mole fraction for Fe in metal

dIW_E_simple = 2*log10(FeO_Mol_E./Fe_Mol_core);
dIW_E_simple(isinf(dIW_E_simple)) = 0;
dIW_E_long = 2*log10((1.148*FeO_Mol_E + 1.319*FeO_Mol_E.^2)./Fe_Mol_core);     %with Rubie+11 param to X for wustite
dIW_E_long(isinf(dIW_E_long)) = 0;

figure(2);
subplot(2,2,3)
hold on
box on
%plot(Accr_model, dIW_E_simple, 'bo', "MarkerSize", 8, "DisplayName", "simple model")
plot(Accr_model, dIW_E_long, 'r.', "MarkerSize", 8, "DisplayName", "Earth calc'd R11 params")
plot(F_temp, dIW_R11, 'k-', "MarkerSize", 8, "DisplayName", "from R11 data")
xlabel("Mass fraction")
ylabel("\DeltaIW")
title("\DeltaIW during accretion")
legend("Location", "southeast")


% Earth core/mantle evolution
if model == 5
    Accr_model = Accr_model*2;
    M_c_ratio = interp1(F_temp, M_c_ratio_temp, Accr_model(1:preGI_idx), 'pchip')/100;
    M_c_ratio_GI = interp1(F_temp, M_c_ratio_temp, 0.99, 'pchip')/100;
    temp_postGI = zeros(1,length(Accr_model)-preGI_idx)+M_c_ratio_GI;
    M_c_ratio = cat(2,M_c_ratio,temp_postGI);
    Accr_model = Accr_model/2;
else
    M_c_ratio = interp1(F_temp, M_c_ratio_temp, Accr_model, 'pchip')/100;
end
M_c = M_E .* M_c_ratio;          %core
M_m = M_E - M_c;                %mantle

FeO_tot = M_m.*OxiWts(5,:)/100; %[kg] mass of FeO for Earth
FeO_added = diff(FeO_tot);      %[kg] difference of FeO mass between impacts

figure(3)
hold on
box on
plot(t, M_E, 'DisplayName', 'Total Mass')
plot(t, M_m, 'DisplayName', 'Mantle Mass')
plot(t, M_c, 'DisplayName', 'Core Mass')
plot(t, FeO_tot, 'DisplayName', 'FeO Mass in Mantle')
ylabel("Mass (kg)")
xlabel("Time (Myr)")
xlim([0 60])
legend('Location', 'northwest')
title("Earth Evolution")

data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpEarly_mantle');
OxiWts_early = data(:,2);
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpLate_mantle');
OxiWts_late = data(:,2);

data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpEarly_core');
MetalWts_early = data(1:6,2);
%M_c_imp_early_ratio = data(7,2)/100;
data = readmatrix('\db\MoleWeights.xlsx', 'Sheet', 'ImpLate_core');
MetalWts_late = data(1:6,2);
%M_c_imp_late_ratio = data(7,2)/100;

M_c_imp = diff(M_c);
M_c_imp_ratio = M_c_imp./M_imp;
M_m_imp = M_imp - M_c_imp;

figure(2)
subplot(2,2,1)
hold on
box on
plot(Accr_model, M_c_ratio*100, 'r.', 'MarkerSize', 8, 'DisplayName', "Earth calc'd")
plot(F_temp, M_c_ratio_temp, 'k-', 'DisplayName', 'From R11 data')
plot(Accr_model(2:end), M_c_imp_ratio*100, 'b.', 'MarkerSize', 8, 'DisplayName', "Imp calc'd")
ylabel("Core Mass (%)")
xlabel("Mass Fraction")
title("Core % Mass vs Accretion")
legend('Location', 'southeast')

FeO_imp = FeO_added./M_m_imp*100;             %[wt %] FeO content of the impactor mantle
FeO_imp_forFig = FeO_imp;
FeO_imp(isnan(FeO_imp)) = 0;
FeO_test = (FeO_imp.*M_m_imp + OxiWts(5,1:end-1).*M_m(1:end-1))./(M_m(2:end));

figure(2)
subplot(2,2,2)
hold on
box on
plot(t(2:end), FeO_test, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Earth FeO')
plot(t(2:end), FeO_imp_forFig, 'b.', 'MarkerSize', 8, 'DisplayName', 'Impactor FeO')
ylabel('FeO Wt %')
xlabel('Time (Myr)')
xlim([0 60])
title("FeO Wt % in Mantles")
legend('Location', 'southeast')


Si_imp = zeros(1, length(FeO_imp));
FeO_Mol_imp = zeros(1, length(FeO_imp));
OxiWts_imp = zeros(length(OxiWts_early), length(FeO_imp));
FeMol_imp = zeros(1, length(FeO_imp));

tempMol = [];   %reset these
MolSum =[];

if model == 5
    Accr_model = 2*Accr_model;
    for i = 1:length(FeO_imp)
        if FeO_imp(i) > 0 && Accr_model(i+1) <= 1.0
            OxiWts_imp(:,i) = OxiWts_early;
            OxiWts_imp(1,i) = 100 - (OxiWts_early(3)+OxiWts_early(7)+OxiWts_early(9)+FeO_imp(i));  %SiO2 content
            Si_imp(i) = OxiWts_imp(1,i);
            OxiWts_imp(5,i) = FeO_imp(i);                                                                %FeO wt% content
            tempMol = OxiWts_imp(:,i)./OxiMolWts;
            MolSum = sum(tempMol);    %total moles
            Mol = tempMol/MolSum;     %mole fraction
            FeO_Mol_imp(i) = Mol(5);      %extract FeO mole fraction

            tempMol = MetalWts_early./MetalMolWts;
            MolSum = sum(tempMol,1);
            Mol_imp_core = tempMol./MolSum;
            FeMol_imp(i) = Mol_imp_core(1);

        elseif FeO_imp(i) > 0 && Accr_model(i+1) > 1.0     %GI as similar comp as Earth
            OxiWts_imp(:,i) = OxiWts(:,i+1);     
            OxiWts_imp(5,i) = FeO_imp(i);
            OxiWts_imp(1,i) = 100 - sum(OxiWts(2:end,i+1));  %SiO2 content
            Si_imp(i) = OxiWts_imp(1,i);                     
            tempMol = OxiWts_imp(:,i)./OxiMolWts;
            MolSum = sum(tempMol);
            Mol = tempMol/MolSum;
            FeO_Mol_imp(i) = Mol(5);

            tempMol = MetalWts./MetalMolWts;
            MolSum = sum(tempMol,1);
            Mol_imp_core = tempMol./MolSum;
            FeMol_imp(i) = Mol_imp_core(1);
        end
    end
    Accr_model = 0.5*Accr_model;
else
    for i = 1:length(FeO_imp)
        if FeO_imp(i) > 0 && Accr_model(i+1) < 0.67
            OxiWts_imp(:,i) = OxiWts_early;
            OxiWts_imp(1,i) = 100 - (OxiWts_early(3)+OxiWts_early(7)+OxiWts_early(9)+FeO_imp(i));  %SiO2 content
            Si_imp(i) = OxiWts_imp(1,i);
            OxiWts_imp(5,i) = FeO_imp(i);                                                                %FeO wt% content
            tempMol = OxiWts_imp(:,i)./OxiMolWts;
            MolSum = sum(tempMol);    %total moles
            Mol = tempMol/MolSum;     %mole fraction
            FeO_Mol_imp(i) = Mol(5);      %extract FeO mole fraction

            tempMol = MetalWts_early./MetalMolWts;
            MolSum = sum(tempMol,1);
            Mol_imp_core = tempMol./MolSum;
            FeMol_imp(i) = Mol_imp_core(1);

        elseif FeO_imp(i) > 0 && Accr_model(i+1) < 1.0
            OxiWts_imp(:,i) = OxiWts_late;
            OxiWts_imp(1,i) = 100 - (OxiWts_late(3)+OxiWts_late(7)+OxiWts_late(9)+FeO_imp(i));  %SiO2 content
            Si_imp(i) = OxiWts_imp(1,i);
            OxiWts_imp(5,i) = FeO_imp(i);                                                       %FeO wt% content
            tempMol = OxiWts_imp(:,i)./OxiMolWts;
            MolSum = sum(tempMol);    %total moles
            Mol = tempMol/MolSum;     %mole fraction
            FeO_Mol_imp(i) = Mol(5);      %extract FeO mole fraction

            tempMol = MetalWts_late./MetalMolWts;
            MolSum = sum(tempMol,1);
            Mol_imp_core = tempMol./MolSum;
            FeMol_imp(i) = Mol_imp_core(1);
        end
    end
end

dIW_imp_simple = 2*log10(FeO_Mol_imp./FeMol_imp);
dIW_imp_simple(isinf(dIW_imp_simple)) = 0;
dIW_imp_long = 2*log10((1.148*FeO_Mol_imp + 1.319*FeO_Mol_imp.^2)./FeMol_imp);     %with Rubie+11 param to X for wustite
dIW_imp_long(isinf(dIW_imp_long)) = 0;

figure(2);
subplot(2,2,3)
hold on
box on
%plot(Accr_model(2:end), dIW_imp_simple, 'b.', "MarkerSize", 8, "DisplayName", "simple model")
plot(Accr_model(2:end), dIW_imp_long, 'b.', "MarkerSize", 8, "DisplayName", "Impactor calc'd R11 params")
legend("Location", "southeast")

%mean values from impactor T and P equilibrium of Table S3a, Rubie+2011
T_early = 2226;
P_early = 5.4;
T_late = 2324;
P_late = 8.8;
T_Elike = 2834;
P_Elike = 42.1;

PV_table = readmatrix('\db\PVcalc_old.xlsx');
P_test = PV_table(:,1);
T_test = [2000, 2500, 3000, 3500, 4000, 4500];
PV_test = PV_table(:,2:end);

PV_temp = interp2(T_test,P_test,PV_test,T_early,P_early);
PV_term_early = PV_temp/log(10);
PV_temp = interp2(T_test,P_test,PV_test,T_late,P_late);
PV_term_late = PV_temp/log(10);
PV_temp = interp2(T_test,P_test,PV_test,T_Elike,P_Elike);
PV_term_Elike = PV_temp/log(10);

FeRatio_imp = zeros(1,length(FeO_imp));              % determine Fe3+/sumFe of impactor based on P, T, dIW
if model == 5
    Accr_model = 2*Accr_model;
    for i = 1:length(FeO_imp)
        if FeO_imp(i) > 0 && Accr_model(i+1) < 0.67
            FeRatio_imp(i) = calcFeRatio_byIW(T_early, P_early, dIW_imp_long(i), PV_term_early, OxiWts_imp(:,i), OxiMolWts_byM);
        elseif FeO_imp(i) > 0 && Accr_model(i+1) <= 1.0
            FeRatio_imp(i) = calcFeRatio_byIW(T_late, P_late, dIW_imp_long(i), PV_term_late, OxiWts_imp(:,i), OxiMolWts_byM);
        elseif FeO_imp(i) > 0 && Accr_model(i+1) > 1.0
            FeRatio_imp(i) = calcFeRatio_byIW(T_Elike, P_Elike, dIW_imp_long(i), PV_term_Elike, OxiWts_imp(:,i), OxiMolWts_byM);
        end
    end
    Accr_model = 0.5*Accr_model;
else
    for i = 1:length(FeO_imp)
        if FeO_imp(i) > 0 && Accr_model(i+1) < 0.67
            FeRatio_imp(i) = calcFeRatio_byIW(T_early, P_early, dIW_imp_long(i), PV_term_early, OxiWts_imp(:,i), OxiMolWts_byM);
        elseif FeO_imp(i) > 0 && Accr_model(i+1) < 1.0
            FeRatio_imp(i) = calcFeRatio_byIW(T_late, P_late, dIW_imp_long(i), PV_term_late, OxiWts_imp(:,i), OxiMolWts_byM);
        end
    end
end


if write == 1
    writematrix(Accr_model, fileOut, 'Sheet', sheetNameImp)
    writematrix(M_c_imp_ratio, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')
    writematrix(Si_imp, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')
    writematrix(FeO_imp, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')
    writematrix(FeRatio_imp, fileOut, 'Sheet', sheetNameImp, 'WriteMode', 'append')

    writematrix(Accr_model, fileOut, 'Sheet', sheetNameE)
    writematrix(M_c_ratio, fileOut, 'Sheet', sheetNameE, 'WriteMode', 'append')
    writematrix(OxiWts, fileOut, 'Sheet', sheetNameE, 'WriteMode', 'append')
end

FeRatio_imp(FeRatio_imp==0) = NaN;   % so it doesn't plot the zeros

data = readmatrix('\db\Compositions_old.xlsx', 'Sheet', 'H04_imp');
FeRatio_imp_old = data(3,:);
FeO_imp_old = data(2,:);
FeRatio_imp_old(FeRatio_imp_old==0) = NaN;
FeO_imp_old(FeO_imp_old==0) = NaN;

figure(2);
subplot(2,2,4)
hold on
box on
plot(t(2:end), FeRatio_imp, 'b.', "MarkerSize", 8, 'DisplayName', 'New Calc')
plot(t(2:end), FeRatio_imp_old, 'ko', 'MarkerSize', 8, 'DisplayName', 'Old Calc')
ylabel("Fe^{3+}/\SigmaFe")
xlabel("Time (Myr)")
xlim([0 60])
title("Impactor Fe^{3+}/\SigmaFe")
legend('Location', 'southeast')

figure(2);
subplot(2,2,2)
hold on
box on
plot(t(2:end), FeO_imp_old, 'ko', 'MarkerSize', 8, 'DisplayName', 'Imp FeO - Old Calc')


