% EL
% Aug 31, 2023

% Overall function to determine geotherm for magma ocean


clear;

% PARAMETERS
model = 5;                      %accretion models, 1 = W90(high), 2 = W90(low), 3 = H00, 4 = H04, 5 = N21
Tp_type = 'Pmo';                %chooose method to calculate Tp, either 'Pmo' or 'U2Q'
    T0 = 1613;                      %[K] initial Tp if using U2Q method
    h = 0.2;                        %energy contribution factor for U

% CONSTANTS
Cp = 1e3;               %[J/kgK] specific heat  
M_E_0 = 5.97e24;        %[kg] present day mass of Earth
M_c_0 = 1.88e24;        %[kg] present day mass of core
rho_imp = 5000;         %[kg/m^3] approximation based on weighted average (0.68Si + 0.32Fe)


[t,Accr_model] = getAccrModel(model);

%determine Earth, impactor, and core mass
% assume core takes up proportional mass of Earth
% note: during impact "n", earth mass is "n" pre-impact and "n+1" post-impact
M_E = M_E_0 * Accr_model;
M_imp = M_E_0 * diff(Accr_model);
M_c = M_E * M_c_0/M_E_0;

% Rubie+2011, eqn 3
rho_m = (4500-3400)*Accr_model+3400;     %scales with planetary mass, use 3400 for upper mantle density for silicate as initial
rho_c = 2.5 * rho_m;

V_c = M_c./rho_c;                           %volume of core
V_m = (M_E - M_c)./rho_m;                   %volume of mantle
R_E = (3/(4*pi)*(V_m + V_c)).^(1/3);        %radius of Earth
R_c = (3/(4*pi)*(M_c./rho_c)).^(1/3);       %radius of core
R_imp = (3/(4*pi)*M_imp/rho_imp).^(1/3);    %radius of impactor 

%test function, pressure at CMB
P_cmb = zeros(1,length(R_c));
for i = 1:length(R_c)
    P_cmb(i) = get2LayerP(rho_m(i), rho_c(i), R_E(i), R_c(i), R_c(i));
end
P_cmb(isnan(P_cmb))=0;

% choose impactor MO depth
P_mo = P_cmb*d_mo_factor;

%corresponding Tp for whole-mantle melting
Tp_cmb = getMOTp(P_cmb/1e9);


% Calculate Tp using Tp_type
% Pmo = MO pressure aligns with liquidus
% U2Q = grav potential to heating MO
switch Tp_type
    case 'Pmo'
        Tp = zeros(1,length(M_E));
        Tp(1) = T0;
        for i = 2:length(M_E)
            if M_imp(i-1) > 0
                Tp(i) = getMOTp(P_mo(i)/1e9);
            else
                Tp(i) = T0;
            end
        end

    case 'U2Q'
        M_imp_temp = cat(2,M_imp,[0]);      %so it matches length of M_E and R_E
        R_imp_temp = cat(2,R_imp,[0]);
        Tp_temp = calcTforU2Q(T0,Cp,M_E,M_c,M_imp_temp,R_E,R_imp_temp,h);

        Tp = zeros(1,length(M_E));
        Tp(1) = T0;
        for i = 2:length(M_E)
            if Tp_temp(i-1) > Tp_cmb(i)
                Tp(i) = Tp_cmb(i);
            else
                Tp(i) = Tp_temp(i-1);
            end
        end

    otherwise
        disp('Could not calculate Tp based on input Tp_type')
end

P_max = round(P_mo,-9);
Tad = zeros((P_max(end)/1e9*2)+1,length(P_max));      %set up Tad and Pad
Pad = zeros((P_max(end)/1e9*2)+1,length(P_max));

% For the first adiabat
Pad_temp = ((0:1e9:(P_max(1)*2))/2)';           %to go by 0.5 GPa
if size(Pad_temp)==1
    Pad_temp = [0, 0.5e9];         %to prevent errors with some of the accretion models
end
Pad_empty = NaN(size(Pad,1)-size(Pad_temp,1),1);
Pad_temp2 = cat(1, Pad_temp, Pad_empty);
Pad = Pad_temp2;

Tad_temp = getMOAdiabat(Tp(1),Pad_temp);
Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
Tad_temp2 = cat(1, Tad_temp, Tad_empty);
Tad = Tad_temp2;


for i = 2:length(P_max)
    if M_imp(i-1) > 0
        Pad_temp = ((0:1e9:(P_max(i)*2))/2)';           %to go by 0.5 GPa
        Pad_empty = NaN(size(Pad,1)-size(Pad_temp,1),1);
        Pad_temp2 = cat(1, Pad_temp, Pad_empty);
        Pad = cat(2,Pad,Pad_temp2);

        Tad_temp = getMOAdiabat(Tp(i),Pad_temp);
        Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);

    else
        Pad_temp = Pad(:,i-1);
        Pad_empty = NaN(size(Pad,1)-size(Pad_temp,1),1);
        Pad_temp2 = cat(1, Pad_temp, Pad_empty);
        Pad = cat(2,Pad,Pad_temp2);

        Tad_temp = Tad(:,i-1);
        Tad_empty = NaN(size(Tad,1)-size(Tad_temp,1),1);
        Tad_temp2 = cat(1, Tad_temp, Tad_empty);
        Tad = cat(2,Tad,Tad_temp2);
    end
end


if write == 1
    writematrix(Pad/1e9, 'geotherms.xlsx')
    writematrix(Tad, 'geotherms.xlsx', 'Range', 'B1') 
end