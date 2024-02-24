% EL
% Sept 26, 2023
%
% Function to calculate the amount of equlibrated silicate
%
% INPUTS:   M_c_imp     [kg] mass of impactor core/Fe, any dimension
%           Eq_type     [] method to calculate equilibrated silicate around droplet
%                           'cyl' = cylindrical tail
%                           'sph' = spherical shell
%
% OUTPUTS:  M_eq        [kg] mass equilibrated in unit time


function [M_eq] = calcEqSi(M_c_imp, Eq_type)

    rho_Fe = 7800;          %[kg/m^3] iron melt density
    r_d = 0.5e-2;           %[m] radius of droplets (~1cm diameter)
    v_d = 0.5;              %[m/s] percolation velocity
    D = 1e-7;               %[m^2/s] chemical diffusity
    rho_Si = 3750;          %[kg/m^3] silicate melt density

    % estimate time evolution of iron on Earth
    dM_Fe = M_c_imp;                %[kg] iron mass accreted per time interval
    m_d = rho_Fe * 4/3*pi*r_d^3;    %[kg] mass of a single droplet
    N_d = dM_Fe/m_d;                % # of droplets per time interval

    % Using eq distance from Rubie+ 2003 (for OLD 'cylindrical tail' modeling)
    d = 1000;               %[m] eq distance from Rubie+ 2003
    dt = d/v_d;             %[s] time interval to fall distance d

    % Ulvrova+2011, diffusion boundary for "drop" with high Re and low
    % viscosity ratio (for NEW 'spherical shell' modeling)
    a = 0.79;
    Pe = r_d*v_d/D;     %Peclet #
    k = 0.1;            %partition coefficient
    Sh = a*Pe^(0.5);    %Sherwood #
    R_D = 1;            %diffusivity ratio ~ 1  (may be 1-10?)
    h = r_d/(a*Pe^(0.5));                             %diffusion boundary
    tau = Pe/3 * ( k/Sh + 1/(10*R_D)) * r_d/v_d;      %time scale for equilibrium

    switch Eq_type
        case 'cyl'                  %cylindrical tail or ring model
            %mass flow rate J = A*v*rho
            J = 2*pi*rho_Si*D*r_d;          %[kg/s] "droplet tail" mass flow rate, based on A=pi*h^2, h=sqrt(2Dr/v)
            %disp(['Cyl mass flux = ', num2str(J*dt)])
            % amount of silicate mass equilibriated in time interval dt (or tau) from iron "rainfall"
            M_eq = J*dt*N_d;                
            
        case 'sph'                   %spherical shell
            J = 2*pi*r_d*h*v_d*rho_Si/2;    %[kg/s] "droplet ring" mass flow rate, divide by 2 for diffusion gradient
            %disp(['Sph mass flux = ', num2str(J*tau)])
            M_eq = J*tau*N_d;               
            
        otherwise
            disp('Could not calculate mass equilibrated based on type chosen')
            M_eq = NaN;
    end
    
end