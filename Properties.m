function [Propiedades] = Properties(Longitudes,Nodos_por_medio,No_elementos,Propiedades_Generales)
% Esta función crea una matriz con todas las propiedades para cada uno de
% los elementos, esto es útil cuando se ensamblan las matrices globales
% de masa y de rigidez.

% Se define una nueva matriz que contiene el número de nodos para cada uno
% de los elementos en el dominio, se hace la suma acumulada de cada uno de
% los nodos.
Nodos_elemento=Longitudes*0;
for i=1:length(Longitudes)
    if i==1
    Contador=Nodos_por_medio(i,1);
    Nodos_elemento(i,1)=Contador;
    elseif i>1
        Contador=Nodos_por_medio(i,1);
        Nodos_elemento(i,1)=Contador+Nodos_elemento(i-1,1);
    end
end

Nodos_elemento=[0;Nodos_elemento];

% Matriz de propiedades del dominio, esta matriz contiene las propiedades
% de todos los elementos, la matriz propiedades generales solo contiene lo
% valores de las propiedades en el dominio, más no las propiedades por elemento
Propiedades=zeros(No_elementos,3);
for i=1:length(Longitudes)
    a=Nodos_elemento(i,1)+1;
    b=Nodos_elemento(i+1,1)+1;
        Propiedades(a:b,1)=Propiedades_Generales(i,1);
        Propiedades(a:b,2)=Propiedades_Generales(i,2);
        Propiedades(a:b,3)=Propiedades_Generales(i,3);
end

end