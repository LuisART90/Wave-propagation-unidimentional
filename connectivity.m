function [Conectividad] = connectivity(No_elementos)
% Se genera una matriz de conectividad, esta matriz permite ensamblar las
% matrices globales, ya que permite identificar los grados de libertad que
% están conectados.
Conectividad=zeros(No_elementos,2);
for i=1:No_elementos
    Conectividad(i,1)=i;
    Conectividad(i,2)=i+1;
end 
end