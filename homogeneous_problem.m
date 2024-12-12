function [outputArg1,outputArg2] = homogeneous_problem(Propiedades_Generales,Conectividad,dx,No_nodos)
%% PRBLEMA HOMOGÉNEO

Modulo_elasticidad=Propiedades_Generales(2,1);
Densidad=Propiedades_Generales(2,2);
Velocidad_Homogeneo=sqrt(Modulo_elasticidad/Densidad);

% Ensamble de la matriz global de rigidez del problema homogéneo
% Matriz de rigidez local
K_Global_Homogeneo=zeros(No_nodos,No_nodos);
for i=1:No_elementos
    K_local_Homogeneo=(Modulo_elasticidad/dx)*[1 -1;-1 1];
    fila=Conectividad(i,:);
    columna=fila;
    K_Global_Homogeneo(fila,columna)=K_Global_Homogeneo(fila,columna)+K_local_Homogeneo;
end

% Ensamble de la matriz global de masa
M_Global_Homogeneo=zeros(No_nodos,No_nodos);
for i=1:No_elementos
    M_local_Homogeneo=(Densidad*dx)*[1/2 0;0 1/2];
    fila=Conectividad(i,:);
    columna=fila;
    M_Global_Homogeneo(fila,columna)=M_Global_Homogeneo(fila,columna)+M_local_Homogeneo;
end


% ENSAMBLE DE LA MATRIZ DE AMORTIGUAMIENTO (RAYLEIGH)
% Aplicando condiciones de frontera
C_Rayleigh_Homogeneo=(alfa*M_Global_Homogeneo)+(beta*K_Global_Homogeneo);
C_absorbente_Homogeneo=zeros(length(C_Rayleigh_Homogeneo),length(C_Rayleigh_Homogeneo));
C_absorbente_Homogeneo(1,1)=Densidad*Velocidad_Homogeneo;
C_absorbente_Homogeneo(end,end)=Densidad*Velocidad_Homogeneo;
C_Global_Homogeneo=C_Rayleigh_Homogeneo+C_absorbente_Homogeneo;

M_Global_Homogeneo(:,Nodos_fijos)=[];
M_Global_Homogeneo(Nodos_fijos,:)=[];

K_Global_Homogeneo(:,Nodos_fijos)=[];
K_Global_Homogeneo(Nodos_fijos,:)=[];

C_Global_Homogeneo(:,Nodos_fijos)=[];
C_Global_Homogeneo(Nodos_fijos,:)=[];

% Se calcula la matriz de masa inversa
Matriz_inversa_Homogeneo=((M_Global_Homogeneo)+((dt)*C_Global_Homogeneo))\eye(length(M_Global_Homogeneo));

Desplazamiento_Homogeneo=zeros(length(M_Global),s);
Fuente_desp_hom=Fuente;
for i=2:s-1
    if i>2
        Fuente_desp_hom=Fuente_desp_hom*0;
    end
    Desplazamiento_Homogeneo(:,i+1)=Matriz_inversa_Homogeneo*(((dt^2)*(Fuente_desp_hom-(K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i))))+((dt)*C_Global_Homogeneo*Desplazamiento_Homogeneo(:,i))+(M_Global_Homogeneo*((2*Desplazamiento_Homogeneo(:,i))-Desplazamiento_Homogeneo(:,i-1))));
end 

% Se obtiene las velocidades
Velocidades_Homogeneo=zeros(length(M_Global_Homogeneo),s);
Fuente_vel_Homogeneo=Fuente;
for i=3:s-2
    if i>3
        Fuente_vel_Homogeneo=Fuente_vel_Homogeneo*0;
    end
    Velocidades_Homogeneo(:,i)=(C_Global_Homogeneo\eye(length(M_Global_Homogeneo)))*(((Fuente_vel_Homogeneo-M_Global_Homogeneo*(Desplazamiento_Homogeneo(:,i+1)-2*Desplazamiento_Homogeneo(:,i)+Desplazamiento_Homogeneo(:,i-1))./(dt^2))-K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i)));
end

Aceleraciones_Homogeneo=zeros(length(M_Global_Homogeneo),s);
Fuente_ace_Homogeneo=Fuente;
for i=3:s-2
    if i>3
        Fuente_ace_Homogeneo=Fuente_ace_Homogeneo*0;
    end 
    Aceleraciones_Homogeneo(:,i)=(M_Global_Homogeneo\eye(length(M_Global_Homogeneo)))*(Fuente_ace_Homogeneo-(C_Global_Homogeneo*Velocidades_Homogeneo(:,i))-(K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i)));
end 
end