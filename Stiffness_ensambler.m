function [K_Global] = Stiffness_ensambler(No_nodos,No_elementos,Propiedades,dx,Conectividad)
% Esta función ensambla la matriz de rigidez de un dominio unidimensional
% ENSAMBLE DE LA MATRIZ GLOBAL DE RIGIDEZ
K_Global=zeros(No_nodos,No_nodos);  % Matriz de rigidez local
K_Local=[1 -1;-1 1];
for i=1:No_elementos
    K_local=(Propiedades(i,1)/dx)*K_Local;
    fila=Conectividad(i,:);
    columna=fila;
    K_Global(fila,columna)=K_Global(fila,columna)+K_local;
end
end