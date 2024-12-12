function [Propiedades_Generales,Longitud,Velocidad] = General_properties(Modulos_elasticidad,Densidades,Longitudes)

% Esta función acomoda las propiedades elásticas e inerciales introducidas 
% en un arreglo para que su identificación sea simple. Se crea un arreglo matricial 
% llamado Propiedades_generales, la primera columna corresponde a cada uno de los
% módulos de elasticidad de los medios que componen el dominio, la segunda
% columna corresponde a cada una de las velocidades, mientras que la
% tercera columna indica la velocidad de propagación de las ondas en cada
% uno de los medios. 

% Se introducen las propiedades del elásticas e inerciales del dominio, la
% primera columna indica los módulos de elasticidad de los medios, mientras
% que la segunda colunma indica las densidades de cada uno de ellos, y la tercera 
% columna indica la velocidad de onda en cada uno de ellos
Propiedades_Generales=zeros(length(Modulos_elasticidad),3);
Propiedades_Generales(:,1)=Modulos_elasticidad';    %Módulos de elasticidad 
Propiedades_Generales(:,2)=Densidades;     %Densidades del dominio
Propiedades_Generales(:,3)=sqrt(Propiedades_Generales(:,1)./Propiedades_Generales(:,2));    %Velocidades del medio

%Longitud total del dominio
Longitud=sum(Longitudes);

% Velocidad de propagación más pequeña presente en el dominio, se usa para
% generar la malla unidimensional
Velocidad=min(Propiedades_Generales(:,3));


end