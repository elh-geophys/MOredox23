% EL
% Aug 30, 2023

% Function to use Turcotte & Schubert (2002), Eqn 2.73 to get the pressure
% at some radius R from center of planet in a 2-layer planet model
%
% INPUTS:  rho_m    [kg/m^3] density of mantle, scalar
%          rho_c    [kg/m^3] density of core, scalar
%          R_E      [m] radius of Earth/planet, scalar
%          R_c      [m] radius of core, scalar
%          R        [m] query radius from center of planet, can be 1D array
%
% OUTPUTS: P        [Pa] pressure as function of radius, can be 1D array

function [P] = get2LayerP(rho_m, rho_c, R_E, R_c, R)
    
    G = 6.6743e-11;         %grav constant
       
    P = zeros(length(R),1);
    
    for i = 1:length(R)

        %T&S Geodyanm Eqn 2.73, Pressure as a function radius from center with varying size of planet
        if R(i) >= R_c && R(i) <= R_E
            %mantle
            P(i) = 4/3*pi*G*rho_m*R(i)^3*(rho_c-rho_m)*(1/R(i)-1/R_E) + 2/3*pi*G*rho_m^2*(R_E^2-R(i)^2); 
        elseif R(i) > 0 && R(i) < R_c
            %core
            P(i) = 0;      %put in eqn
        elseif R(i) < 0
            P(i) = NaN;
            disp("Error: Input radius less than zero")
        elseif R(i) > R_E
            P(i) = NaN;
            disp("Error: Input radius larger than planet radius")
        else
            P(i) = NaN;
            disp("Error: Could not calculate P with given inputs")
        end
    end

end