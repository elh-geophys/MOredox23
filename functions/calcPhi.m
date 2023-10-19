%Calculate the melt fraction (phi) for MO geotherm
%Code from JK
%Parameterized in Korenaga (2023)

function [phi,dphidT] = calcPhi(T,Ts,T40,Tl)

dphidT=0;
if T>Tl
  phi = 1;
elseif T>Ts && T<Tl
  if T>=T40
    phi = 0.4+0.6*(T-T40)/(Tl-T40);
    dphidT = 0.6/(Tl-T40);
  else
    phi = 0.4*(T-Ts)/(T40-Ts);
    dphidT = 0.4/(T40-Ts);
  end
else % T<=Ts
  phi = 0;
end
