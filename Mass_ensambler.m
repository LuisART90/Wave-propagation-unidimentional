function [M_Global] = Mass_ensambler(No_nodos,No_elementos,Propiedades,dx,Conectividad)
% Esta función ensambla la matriz de masa global
M_Global=zeros(No_nodos,No_nodos);
M_Local=[1/2 0;0 1/2];
for i=1:No_elementos
    M_local=(Propiedades(i,2)*dx)*M_Local;
    fila=Conectividad(i,:);
    columna=fila;
    M_Global(fila,columna)=M_Global(fila,columna)+M_local;
end
end