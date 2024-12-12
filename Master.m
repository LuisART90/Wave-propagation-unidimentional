
close all
clear 
clc

% Se introducen los datos
Frecuencia_simulacion=0.5;  % Hz
s=1500;                     % Pasos de tiempo a simular
Distancia_fuente=19;        % Punto de aplicación de la fuente impulsiva

% Factores de amortiguamiento
alfa=0.001;     % Comúnmente se define un valor igual o menor a 0.05
beta=0.001;

% Propiedades geométricas, elásticas e inerciales del dominio
% Propiedades geométricas
Longitudes=[10;10];         % Longitud de cada uno de los medios que componen el dominio
Modulos_elasticidad=[1,8];  % Módulos de elasticidad de cada uno de los medios que componen el dominio
Densidades=[1,2];           % Densidades de cada uno de los medios que componen el dominio

% Se crea un arreglo para guardar las propiedades en todo el dominio
[Propiedades_Generales,Longitud,Velocidad] = General_properties(Modulos_elasticidad,Densidades,Longitudes);

% Se crea la malla para el dominio unidimensional
[dt,dx,Long_onda,No_elementos,No_nodos,Nodos_por_medio] = Parameter_mesh(Velocidad,Frecuencia_simulacion,Longitud,Longitudes);

% Propiedades por elemento
[Propiedades] = Properties(Longitudes,Nodos_por_medio,No_elementos,Propiedades_Generales);

% Conectividad entre elementos
[Conectividad] = connectivity(No_elementos);

% Matriz de grados de libertad, los nodos que posean un valor igual a cero
% serán nodos fijos, nodos cuyos valores sea igual a 1 estarán libres. Al
% final se obtiene una matriz de grados de libertad totales.
Nodos=ones(No_nodos,1);
% % Suponiendo que uno de los nodos de los extremos está fijo
Nodos(1,1)=0;
% Número de grados de libertad por nodo 
GL_nodo=1;
% Se buscan las celdas de las matriz cuyo valor sea igual de cero, es
% decir los nodos fijos
Nodos_fijos=find(Nodos(:,1)==0);

% Se obtienen las matrices de masa y de rigidez globales y de
% amortiguamiento
[K_Global] = Stiffness_ensambler(No_nodos,No_elementos,Propiedades,dx,Conectividad);
[M_Global] = Mass_ensambler(No_nodos,No_elementos,Propiedades,dx,Conectividad);
[C_Global] = Damping_Absorbing_boundaries(alfa,beta,M_Global,K_Global,Propiedades);

% Aplicando condiciones de frontera para resolver el problema
[K_Global,M_Global,C_Global] = Boundary_conditions(M_Global,K_Global,C_Global,Nodos_fijos);

% Matriz de masa inversa
Matriz_inversa=((M_Global)+((dt)*C_Global))\eye(length(M_Global));

% Fuente impulsiva que genera el movimiento
[Fuente] = Source(Distancia_fuente,dx,Frecuencia_simulacion,Longitud,dt,Nodos_fijos);

% Solución del problema (OPCIONAL)
[Desplazamiento,Velocidades,Aceleraciones] = Solver(dt,Fuente,Matriz_inversa,K_Global,M_Global,C_Global,s);

%% Se resuelve el problema homogéneo

