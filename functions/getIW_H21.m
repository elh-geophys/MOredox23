% EL
% Feb 2023
% Updated 09-06-2023
%
% IW Buffer
% Hirschmann 2021
%
% Inputs:   P       [Pa] Pressure 
%           T       [K] temperature
% 
% Output:   IW at each given P-T value as a log10 value
%
% Note, H21 does not recommend using this function directly past 100GPa.
% We've compromised by linearly extrapolating past 100GPa

function [IW] = getIW_H21(P,T)
    
    P = P/1e9;              %H21 did this in GPa, so convert

    %Table 1 from Hirschmann 2021
    m_cc = [6.844864,       1.175691e-1,    1.143873e-3,    0,              0;
            5.791364e-4,    -2.891434e-4,   -2.737171e-7,   0,              0;
            -7.971469e-5,   3.198005e-5,    0,              1.059554e-10,   2.014461e-7;
            -2.769002e4,    5.285977e2,     -2.919275,      0,              0];

    m_hcp = [8.463095,      -3.000307e-3,   7.213445e-5,    0,              0;
             1.148738e-3,   -9.352312e-5,   5.161592e-7,    0,              0;
            -7.448624e-4,   -6.329325e-6,   0,              -1.407339e-10,  1.830014e-4;
            -2.782082e4,    5.285977e2,     -8.473231e-1,   0,              0];

    %for distinction between bcc/fcc iron and hcp,
    %Table 1 note
    x0 = -18.64;
    x1 = 0.04359;
    x2 = -5.069e-6;

    IW = zeros(length(P),length(T));

    for j = 1:length(T)         %T loop
        P_test = x0 + x1*T(j) + x2.*T(j)^2;

        for i = 1:length(P)     %P loop
            if P(i)>P_test           %use hcp formulation
                e = m_hcp(1,1) + m_hcp(1,2)*P(i) + m_hcp(1,3)*P(i)^2 + m_hcp(1,4)*P(i)^3 + m_hcp(1,5)*P(i)^0.5;
                f = m_hcp(2,1) + m_hcp(2,2)*P(i) + m_hcp(2,3)*P(i)^2 + m_hcp(2,4)*P(i)^3 + m_hcp(2,5)*P(i)^0.5;
                g = m_hcp(3,1) + m_hcp(3,2)*P(i) + m_hcp(3,3)*P(i)^2 + m_hcp(3,4)*P(i)^3 + m_hcp(3,5)*P(i)^0.5;
                h = m_hcp(4,1) + m_hcp(4,2)*P(i) + m_hcp(4,3)*P(i)^2 + m_hcp(4,4)*P(i)^3 + m_hcp(4,5)*P(i)^0.5;

                IW(i,j) = e + f*T(j) + g*T(j)*log(T(j)) + h/T(j);

            else                        %use fcc/bcc formulation
                a = m_cc(1,1) + m_cc(1,2)*P(i) + m_cc(1,3)*P(i)^2 + m_cc(1,4)*P(i)^3 + m_cc(1,5)*P(i)^0.5;
                b = m_cc(2,1) + m_cc(2,2)*P(i) + m_cc(2,3)*P(i)^2 + m_cc(2,4)*P(i)^3 + m_cc(2,5)*P(i)^0.5;
                c = m_cc(3,1) + m_cc(3,2)*P(i) + m_cc(3,3)*P(i)^2 + m_cc(3,4)*P(i)^3 + m_cc(3,5)*P(i)^0.5;
                d = m_cc(4,1) + m_cc(4,2)*P(i) + m_cc(4,3)*P(i)^2 + m_cc(4,4)*P(i)^3 + m_cc(4,5)*P(i)^0.5;

                IW(i,j) = a + b*T(j) + c*T(j)*log(T(j)) + d/T(j);
            end
        end
    end

end