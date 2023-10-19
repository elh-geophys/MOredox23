% EL
% RT instability

clear;
clf;

% data from JK adiabat code for PREM
data = readmatrix('geotherms.xlsx', 'Sheet', 'PREM');
P = data(:,2);          %[GPa]
z = data(:,3);          %[km] depth
rho_s = data(:,4);      %[kg/m^3] PREM density

mu_s = 1e15;            
%[Pa s] viscosity of solid underlying layer (somewhere around 10^13 to 10^18)
                        
mu_l = 0.01;             
%[Pa s] viscosity of metal-rich layer
% Assume similar to outer core (deWijs+ 1998)

rho_l = 7800;           
%[kg/m^3] metal-rich layer density (Rubie+ 2003)

h = logspace(0, 5, 6);      
%[m] metal-rich layer thickness
% ~5km stated in Wood+(2006), based on Rubie+(2003)
% Karato&Murthy(1997) estimated max thickness based on accretion & RT to be between 1m to
% 100km (viscosity of 1e15 to 1e23 for bottom layer, respectively)

[~,idx] = min(abs(rho_s-rho_l));            %to determine where the density contrast changes

% ------with Mondal&Korenaga growth rate
a0 = -10.13;     %the 10 term
b1 = 1.86;      %the R1 term
c1 = 1.3;       %the delta-rho term
c2 = 0.87;      %the rho1 term
b2 = 0.047;     %the delta-R term
d1 = -0.98;     %the mu1 term
d2 = -0.054;    %the mu2 term

R1 = (6371-z)*1e3;
R2 = R1+h;
rho1 = rho_s;
rho2 = rho_l;
mu1 = mu_s;
mu2 = mu_l;

sigma = 10^a0 .* R1.^b1 .* (rho2-rho1).^c1 .* rho1.^c2 .* (R2-R1).^b2 .* mu1.^d1 .* mu2.^d2;

tau = sigma.^-1;
tau_dy = tau/(60*60*24);
tau_hr = tau/(60*60);

newcolors = 1/255*[0 0 0; ...
                   40 40 40; ...
                   80 80 80; ...
                   120 120 120;
                   160 160 160;
                   200 200 200;
                   240 240 240];

figure(1);
subplot(1, 2, 1)
colororder(newcolors)
plot(P(1:idx,:), sigma(1:idx,:), 'LineWidth', 1.5)
xlabel("Pressure [GPa]")
ylabel("Growth rate, \sigma [s^{-1}]")
xlim([min(P), max(P)])
legend(["10^0 m", "10^1 m", "10^2 m", "10^3 m", "10^4 m", "10^5 m"], "Location", "northeast")

subplot(1, 2, 2)
colororder(newcolors)
plot(P(1:idx,:), tau_hr(1:idx,:), 'LineWidth', 1.5)
xlabel("Pressure [GPa]")
ylabel("Instability time, \sigma^{-1} [hours]")
xlim([min(P), max(P)])
ylim([0, 30])
