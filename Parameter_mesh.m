function [dt,dx,Long_onda,No_elementos,No_nodos,Nodos_por_medio] = Parameter_mesh(Velocidad,Frecuencia_simulacion,Longitud,Longitudes)


    % En esta función se definen los parámetros para construir la malla de
    % elementos finitos 
    %Se definen los parámetros para definir la malla del dominio
    No_ptos_long_onda=16;   % Para emplear elementos de interpolación lineales
    Long_onda=Velocidad/Frecuencia_simulacion;
    dx=Long_onda/No_ptos_long_onda;
    No_elementos=round(Longitud/dx);
    No_nodos=No_elementos+1;
    Nodos_por_medio=round(Longitudes./dx);
    
    C=0.5;     % Condición de estabilidad
    dt=C*dx/Velocidad;

end