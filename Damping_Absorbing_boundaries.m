function [C_Global] = Damping_Absorbing_boundaries(alfa,beta,M_Global,K_Global,Propiedades)
% Ensamble de la matriz de amortoguamiento de Rayleigh
C_Rayleigh=(alfa*M_Global)+(beta*K_Global);
% Fronteras absorbentes
C_absorbente=zeros(length(C_Rayleigh),length(C_Rayleigh));
C_absorbente(1,1)=Propiedades(1,2)*Propiedades(1,3);
C_absorbente(end,end)=Propiedades(end,2)*Propiedades(end,3);
% Matriz de amortiguamiento
C_Global=C_Rayleigh+C_absorbente;
end