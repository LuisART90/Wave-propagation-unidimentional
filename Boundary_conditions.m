function [K_Global,M_Global,C_Global] = Boundary_conditions(M_Global,K_Global,C_Global,Nodos_fijos)
% Esta función identifica los grados de libertad del dominio y elimina los
% grados de libertad fijos de las matrices de masa, de rigidez y de
% amortiguamiento

% % Aplicando condiciones de frontera
M_Global(:,Nodos_fijos)=[];
M_Global(Nodos_fijos,:)=[];

K_Global(:,Nodos_fijos)=[];
K_Global(Nodos_fijos,:)=[];

C_Global(:,Nodos_fijos)=[];
C_Global(Nodos_fijos,:)=[];

end