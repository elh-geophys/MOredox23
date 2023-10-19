% EL
% Aug 30, 2023
%
% Function to get adiabat in magma ocean with a given potential
% temperature. Based on parameterization from Korenaga (2023), which is
% from Miyazaki & Korenaga (2019b).
%
% NOTE:  Only gives adiabat for above liquidus. For more extensive
% geotherm, use JK code.
%
% INPUTS:   Tp      [K] potential temperature
%           P       [Pa] Pressure profile, 1D array
%
% OUTPUTS:  Tad     [K] adiabat profile as function of P, 1D array
%

function [Tad] = getMOAdiabat(Tp, P, data)

    %data = readmatrix('\db\geotherms_combo.xlsx');
    P_data = data(:,1)*1e9;
    Tsol_data = data(:,3);
    Tliq_data = data(:,4);
    T40_data = data(:,5);
    
    Tsol = interp1(P_data, Tsol_data, P);
    Tliq = interp1(P_data, Tliq_data, P);
    T40 = interp1(P_data, T40_data, P);

    Tad = zeros(length(P),1);
    Tad(1) = Tp;
    
%     figure(3);
%     hold on
%     box on
%     plot(Tsol,P/1e9, 'k-')
%     plot(Tliq,P/1e9, 'k-')
%     plot(T40,P/1e9, 'k:')
    
    for j=2:length(P)
        dp = P(j)-P(j-1);
        
        [phi1,dphidT] = calcPhi(Tad(j-1),Tsol(j-1),T40(j-1),Tliq(j-1));
        Hf = interp1([0 max(P)],[6e5 9e6],P(j));
  
        alpha = 3.622e-5*exp(-2.377e-5*Tad(j-1)-0.0106*P(j-1)/1e9);
        cp =  627+0.441*Tad(j-1)-0.211*P(j-1)/1e9;
        rho = 2870-0.082*Tad(j-1)+162*(P(j-1)/1e9)^0.58;
        dTdp = alpha*Tad(j-1)/(rho*cp);      
        dT = dTdp*dp;
        Tad0 = Tad(j-1)+dT;
        
        if Tad0<Tliq(j) && Tad0>Tsol(j)    %bewteen liquidus & solidus                        %between solidus & liquidus
            [~,dphidT1] = calcPhi(0.5*(T40(j)+Tliq(j)),...
                         Tsol(j),T40(j),Tliq(j));
            [~,dphidT2] = calcPhi(0.5*(T40(j)+Tsol(j)),...
                         Tsol(j),T40(j),Tliq(j));
            if phi1>0.4                                         
                phi2 = (Hf/cp*phi1+1/dphidT1-Tliq(j)+Tad0)/(Hf/cp+1/dphidT1);
                
                if phi2>0.4
                    newT = Tliq(j)-1/dphidT1*(1-phi2);
                else
                    newT = T40(j)-1/dphidT2*(0.4-phi2);
                end
                
            else            %no calc beyond phi < 0.4
              %phi2 = (Hf/cp*phi1+1/dphidT2*0.4-T40(j)+Tad0)...
              %   /(Hf/cp+1/dphidT2);
              %newT = T40(j)-1/dphidT2*(0.4-phi2);
              newT = NaN;
              
            end

            Tad(j) = newT;

        else                           %else standard adiabat
            Tad(j) = Tad0;
        
        end
        
    end
%     figure(3);
%     plot(Tad,P/1e9,'b', "LineWidth",1.5)
% 
%     pause(0.3)
end
