% EL
% Earth accretion models
% December 2022
%
% type = 1 is Wetherill (1990) with growth factor 10 Myr
% type = 2 is Wetherill (1990) with growth factor 13 Myr
% type = 3 is Halliday (2000) with large GI at 50 Myr
% type = 4 is Halliday (2004) with discrete impacts and GI at 55 Myr
% type = 5 is Nesvorny (2021) with GI at 41 Myr with 1:1 mass ratio


function [t, model] = getAccrModel(type)

    t = linspace(0, 100, 1000);          % [Myr], where 0 = Earth origin

    switch type
        case 1
            % Wetherill (1990) with growth factor 10 Myr
            lambda_W90_high = 1/10;
            f_W90_high = 1 - exp(-lambda_W90_high*t);

            model = f_W90_high;

        case 2
            % Wetherill (1990) with growth factor 13 Myr
            lambda_W90_mid = 1/13;
            f_W90_mid = 1 - exp(-lambda_W90_mid*t);

            model = f_W90_mid;

        case 3
            % Halliday (2000) with large GI at 50 Myr
            t_GI = 50;                  % [Myr] approx time of giant impact (GI)
            f_postGI = 0.90;            % frac of earth mass post-GI (Earth + Theia mass)
            f_GI = 7/10*f_postGI;       % frac of earth mass pre-GI (to keep 7:3 ratio)
            f_IMP = f_postGI - f_GI;    % frac of mass added during GI
            lambda1 = (1/t_GI * log(1/(1-f_GI)));
            lambda2 = (1/t_GI * log(1/(1-(f_GI + f_IMP))));

            f_H00 = zeros(1, length(t));
            for i = 1:length(t)
               if t(i) <= t_GI
                   f_H00(i) = 1 - exp(-lambda1*t(i));
               else
                   f_H00(i) = 1 - exp(-lambda2*t(i));
               end
            end

            model = f_H00;
        
        case 4
            % Halliday (2004) with discrete impacts and GI at 55 Myr
            % this version follows Wetherill(1990) with growth factor 19 Myr
            lambda_W90_low = 1/19;
            f_W90_low = 1 - exp(-lambda_W90_low*t);

            f_H04_low = zeros(1, length(t));
            f_H04_low(1) = 0.01;
            for i = 2:length(t)
                if round(f_H04_low(i-1),3) < 0.1
                    % 1% increase 
                    [~,idx] = min(abs(f_W90_low-(f_H04_low(i-1)+0.005)));
                    if i < idx
                        f_H04_low(i) = f_H04_low(i-1);
                    else
                        f_H04_low(i) = f_H04_low(i-1) + 0.01;
                    end
                elseif round(f_H04_low(i-1),3) <0.3
                    % 2% increase
                    [~,idx] = min(abs(f_W90_low-(f_H04_low(i-1)+0.01)));
                    if i < idx
                        f_H04_low(i) = f_H04_low(i-1);
                    else
                        f_H04_low(i) = f_H04_low(i-1) + 0.02;
                    end
                elseif round(f_H04_low(i-1),3) < 0.9
                    % 4% increase
                    [~,idx] = min(abs(f_W90_low-(f_H04_low(i-1)+0.02)));
                    if i < idx
                        f_H04_low(i) = f_H04_low(i-1);
                    else
                        f_H04_low(i) = f_H04_low(i-1) + 0.04;
                    end 
                else
                    % 9% increase for GI
                    [~,idx] = min(abs(f_W90_low-(f_H04_low(i-1)+0.045)));
                    if i <= idx
                        f_H04_low(i) = f_H04_low(i-1);
                    else
                        f_H04_low(i) = f_H04_low(i-1) + 0.09;
                    end 
                end
            end

            model = f_H04_low;

        case 5
            % Nesvorny+ 2021, fast initial growth, large GI at 41 Myr
            t_GI = 41;                  % [Myr] approx time of giant impact (GI)

            f_N21 = zeros(1, length(t));
            for i = 1:length(t)
               if t(i) <= 1
                   f_N21(i) = 0.1;
               elseif t(i) <= 2
                   f_N21(i) = 0.2;
               elseif t(i) <= 3
                   f_N21(i) = 0.3;
               elseif t(i) <= 4
                   f_N21(i) = 0.4;
               elseif t(i) <= t_GI
                   f_N21(i) = 0.495;
               else
                   f_N21(i) = 0.99;
               end
            end

            model = f_N21;
        
        otherwise
            model = [];
            disp("Error: input 'type' was not recognized as a corresponding model");
    end

end