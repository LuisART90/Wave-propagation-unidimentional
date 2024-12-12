

%Programa de elementos finitos unidimensionales, se aplica el método de
%reduccíón del dominio
close all
clear 
clc

% profile on
% El medio geológico original se define como heterogéneo, tres medios
% distintos
Frecuencia_simulacion=1;    %Hz
s=3000;     %Pasos de tiempo a simular
Longitudes=[10;10];
Longitud=sum(Longitudes);
Distancia_fuente=19;

% Factores de amortiguamiento
alfa=0.005;
beta=0.005;

% Se introducen las propiedades del elásticas e inerciales del dominio, la
% primera columna indica los módulos de elasticidad de los medios, mientras
% que la segunda colunma indica las densidades de cada uno de ellos, y la tercera 
% columna indica la velocidad de onda en cada uno de ellos
Propiedades_Generales=zeros(length(Longitudes),3);
Propiedades_Generales(:,1)=[1;8];    %Módulos de elasticidad 
Propiedades_Generales(:,2)=[1;2];     %Densidades del dominio
Propiedades_Generales(:,3)=sqrt(Propiedades_Generales(:,1)./Propiedades_Generales(:,2));    %Velocidades del medio
Velocidad=Propiedades_Generales(2,3);

%Se definen los parámetros para definir la malla del dominio
No_ptos_long_onda=16;   % Para emplear elementos de interpolación lineales
Long_onda=Velocidad/Frecuencia_simulacion;
dx=Long_onda/No_ptos_long_onda;
No_elementos=round(Longitud/dx);
No_nodos=No_elementos+1;    % lineales
% No_nodos=(2*No_elementos)+1;    % cuadráticos
Nodos_por_medio=round(Longitudes./dx);

C=0.25;     % Condición de estabilidad
dt=C*dx/Velocidad;

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

% Se genera una matriz de conectividad con elementos lineales
Conectividad=zeros(No_elementos,2);
for i=1:No_elementos
    Conectividad(i,1)=i;
    Conectividad(i,2)=i+1;
end 

% % Matriz de conectividad con elementos cuadráticos
% Conectividad=zeros(No_elementos,3);
% for i=1:No_elementos
%     Conectividad(i,1)=(2*i)-1;
%     Conectividad(i,2)=(2*i);
%     Conectividad(i,3)=(2*i)+1;
% end 

% Matriz de grados de libertad, los nodos que posean un valor igual a cero
% serán nodos fijos, nodos cuyos valores sea igual a 1 estarán libres. Al
% final se obtiene una matriz de grados de libertad totales.
Nodos=ones(No_nodos,1);
% Suponiendo que los nodos de los extremos están fijos
Nodos(1,1)=0;
Nodos(end,1)=0;

% Número de grados de libertad por nodo 
GL_nodo=1;

% Se buscan las celdas de las matriz cuyo valor sea igual de cero, es
% decir los nodos fijos
Nodos_fijos=find(Nodos(:,1)==0);

% ENSAMBLE DE LA MATRIZ GLOBAL DE RIGIDEZ
K_Global=zeros(No_nodos,No_nodos);  % Matriz de rigidez local
K_Local=[1 -1;-1 1];    % lineales
% K_Local=(1/3)*[7 -8 1;-8 16 -8;1 -8 7]; % Cuadráticos
for i=1:No_elementos
    K_local=(Propiedades(i,1)/dx)*K_Local;
    fila=Conectividad(i,:);
    columna=fila;
    K_Global(fila,columna)=K_Global(fila,columna)+K_local;
end

% EMSAMBLE DE LA MATRIZ DE GLOBAL DE MASA
M_Global=zeros(No_nodos,No_nodos);
M_Local=[1/2 0;0 1/2];  % lineales
% M_Local=(1/30)*[5 0 0;0 20 0;0 0 5];   % Cuadráticos
for i=1:No_elementos
    M_local=(Propiedades(i,2)*dx)*M_Local;
    fila=Conectividad(i,:);
    columna=fila;
    M_Global(fila,columna)=M_Global(fila,columna)+M_local;
end

% ENSAMBLE DE LA MATRIZ DE AMORTIGUAMIENTO (RAYLEIGH)
C_Rayleigh=(alfa*M_Global)+(beta*K_Global);
C_absorbente=zeros(length(C_Rayleigh),length(C_Rayleigh));
C_absorbente(2,2)=Propiedades(1,2)*Propiedades(1,3);
C_absorbente(end-1,end-1)=Propiedades(end,2)*Propiedades(end,3);
C_Global=C_Rayleigh+C_absorbente;

% % Aplicando condiciones de frontera
M_Global(:,Nodos_fijos)=[];
M_Global(Nodos_fijos,:)=[];

K_Global(:,Nodos_fijos)=[];
K_Global(Nodos_fijos,:)=[];

C_Global(:,Nodos_fijos)=[];
C_Global(Nodos_fijos,:)=[];

Matriz_inversa=((M_Global)+((dt/2)*C_Global))\eye(length(M_Global));

% Función Gaussiana 
mu=Distancia_fuente;
sigma=dx;
if Frecuencia_simulacion>1
x=0:dx:Longitud;
elseif Frecuencia_simulacion<=1
    x=0:dx:Longitud;
end
Gaussian=((1/(sigma*sqrt(2*pi)))*exp(-(((x-mu).^2)/(2*sigma*sigma))))';
Fuente=Gaussian*10;
Fuente(Nodos_fijos,:)=[];

% Fuente=zeros(No_nodos-2,1);
% Fuente(58,1)=50;
 
% Desplazamiento=zeros(length(M_Global),s);
% Fuente_desp=Fuente;
% for i=2:s-1
%     if i>2
%         Fuente_desp=Fuente*0;
%     end
%     Desplazamiento(:,i+1)=Matriz_inversa*(((dt^2)*(Fuente_desp-(K_Global*Desplazamiento(:,i))))+((dt/2)*C_Global*Desplazamiento(:,i-1))+(M_Global*((2*Desplazamiento(:,i))-Desplazamiento(:,i-1))));
% end 
% 
% % Se obtiene las velocidades
% Velocidades=zeros(length(M_Global),s);
% Fuente_vel=Fuente;
% for i=2:s-1
%     if i>2
%         Fuente_vel=Fuente_vel*0;
%     end
%     Velocidades(:,i)=(C_Global\eye(length(M_Global)))*(((Fuente_vel-M_Global*(Desplazamiento(:,i+1)-2*Desplazamiento(:,i)+Desplazamiento(:,i-1))./(dt^2))-K_Global*Desplazamiento(:,i)));
% end
% 
% Aceleraciones=zeros(length(M_Global),s);
% Fuente_ace=Fuente;
% for i=1:s
%     if i>1
%         Fuente_ace=Fuente_ace*0;
%     end 
%     Aceleraciones(:,i)=(M_Global\eye(length(M_Global)))*(Fuente_ace-(C_Global*Velocidades(:,i))-(K_Global*Desplazamiento(:,i)));
% end 
% 
% % VISUALIZACIÓN DE LOS DESPLAZAMIENTOS, VELOCIDAD Y ACELERACIÓN
% close all
% figure
% for i=1:s
% 
% subplot(3,1,1)
% plot(Desplazamiento(:,i))
% axis([-10 length(M_Global)+10 min(min(Desplazamiento)) max(max(Desplazamiento))])
% title('Desplazamiento')
% grid on
% 
% subplot(3,1,2)
% plot(Velocidades(:,i))
% axis([-10 length(M_Global)+10 min(min(Velocidades)) max(max(Velocidades))])
% title('Velocidad')
% grid on
% 
% subplot(3,1,3)
% plot(Aceleraciones(:,i))
% axis([-10 length(M_Global)+10 min(min(Aceleraciones)) max(max(Aceleraciones))])
% title('Aceleración')
% grid on
% 
% drawnow
% end

% % OBTENCIÓN DE LOS VALORES PARA DEFINIR LA MATRIZ DE AMORTIGUAMIENTO
% % Valores y vectores propios
% % los valores propios son las frecuencias naturales del sistema asociada a
% % cada uno de los modos de vibrar
% [Modos_de_vibrar,frecuencias_naturales]=eig(K_Global);
% 
% % Se visualizan los modos de vibrar del sistema
% % figure
% % for i=1:length(K_Global)
% %     plot(Modos_de_vibrar(:,i))
% %     title('Modos de vibrar')
% %     axis([-10 No_nodos+10 -1 1])
% %     grid on
% %     drawnow
% %     pause(1)
% % end 
% 
% % Coeficiente de amortiguamiento global
% Zeta=0.05;
% 
% Modo_1=1;
% Modo_2=length(M_Global);
% 
% frecuencia_1=frecuencias_naturales(Modo_1,Modo_1);
% frecuencia_2=frecuencias_naturales(Modo_2,Modo_2);
% 
% alfa=(2*Zeta*frecuencia_1*frecuencia_2)/(frecuencia_1+frecuencia_2);
% beta=(2*Zeta)/(frecuencia_1+frecuencia_2);
% 
% figure
% for i=1:s
%     plot(Modos_de_vibrar(:,i))
%     drawnow
%     pause
% end 
% 



%% PRBLEMA HOMOGÉNEO

Modulo_elasticidad=Propiedades_Generales(2,1);
Densidad=Propiedades_Generales(2,2);
Velocidad_Homogeneo=sqrt(Modulo_elasticidad/Densidad);

% Ensamble de la matriz global de rigidez del problema homogéneo
% Matriz de rigidez local
K_Global_Homogeneo=zeros(No_nodos,No_nodos);
for i=1:No_elementos
    K_local_Homogeneo=(Modulo_elasticidad/dx)*[1 -1;-1 1];
%     K_local_Homogeneo=(Modulo_elasticidad/dx)*(1/3)*[7 -8 1;-8 16 -8;1 -8 7];   % cuadráticos
    fila=Conectividad(i,:);
    columna=fila;
    K_Global_Homogeneo(fila,columna)=K_Global_Homogeneo(fila,columna)+K_local_Homogeneo;
end

% Ensamble de la matriz global de masa
M_Global_Homogeneo=zeros(No_nodos,No_nodos);
for i=1:No_elementos
    M_local_Homogeneo=(Densidad*dx)*[1/2 0;0 1/2];
%     M_local_Homogeneo=(Densidad*dx)*(1/30)*[5 0 0;0 20 0;0 0 5];   % cuadráticos
    fila=Conectividad(i,:);
    columna=fila;
    M_Global_Homogeneo(fila,columna)=M_Global_Homogeneo(fila,columna)+M_local_Homogeneo;
end


% ENSAMBLE DE LA MATRIZ DE AMORTIGUAMIENTO (RAYLEIGH)
% Aplicando condiciones de frontera
C_Rayleigh_Homogeneo=(alfa*M_Global_Homogeneo)+(beta*K_Global_Homogeneo);
C_absorbente_Homogeneo=zeros(length(C_Rayleigh_Homogeneo),length(C_Rayleigh_Homogeneo));
C_absorbente_Homogeneo(2,2)=Densidad*Velocidad_Homogeneo;
C_absorbente_Homogeneo(end-1,end-1)=Densidad*Velocidad_Homogeneo;
C_Global_Homogeneo=C_Rayleigh_Homogeneo+C_absorbente_Homogeneo;

M_Global_Homogeneo(:,Nodos_fijos)=[];
M_Global_Homogeneo(Nodos_fijos,:)=[];

K_Global_Homogeneo(:,Nodos_fijos)=[];
K_Global_Homogeneo(Nodos_fijos,:)=[];

C_Global_Homogeneo(:,Nodos_fijos)=[];
C_Global_Homogeneo(Nodos_fijos,:)=[];

% Se calcula la matriz de masa inversa
Matriz_inversa_Homogeneo=((M_Global_Homogeneo)+((dt/2)*C_Global_Homogeneo))\eye(length(M_Global_Homogeneo));

Desplazamiento_Homogeneo=zeros(length(M_Global),s);
Fuente_desp_hom=Fuente;
for i=2:s-1
    if i>2
        Fuente_desp_hom=Fuente_desp_hom*0;
    end
    Desplazamiento_Homogeneo(:,i+1)=Matriz_inversa_Homogeneo*(((dt^2)*(Fuente_desp_hom-(K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i))))+((dt/2)*C_Global_Homogeneo*Desplazamiento_Homogeneo(:,i-1))+(M_Global_Homogeneo*((2*Desplazamiento_Homogeneo(:,i))-Desplazamiento_Homogeneo(:,i-1))));
end 

% Se obtiene las velocidades
Velocidades_Homogeneo=zeros(length(M_Global_Homogeneo),s);
Fuente_vel_Homogeneo=Fuente;
for i=2:s-1
    if i>2
        Fuente_vel_Homogeneo=Fuente_vel_Homogeneo*0;
    end
    Velocidades_Homogeneo(:,i)=(C_Global_Homogeneo\eye(length(M_Global_Homogeneo)))*((Fuente_vel_Homogeneo-M_Global_Homogeneo*(Desplazamiento_Homogeneo(:,i+1)-2*Desplazamiento_Homogeneo(:,i)+Desplazamiento_Homogeneo(:,i-1))/dt^2)-(K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i)));
end

Aceleraciones_Homogeneo=zeros(length(M_Global_Homogeneo),s);
Fuente_ace_Homogeneo=Fuente;
for i=1:s
    if i>1
        Fuente_ace_Homogeneo=Fuente_ace_Homogeneo*0;
    end 
    Aceleraciones_Homogeneo(:,i)=(M_Global_Homogeneo\eye(length(M_Global_Homogeneo)))*(Fuente_ace_Homogeneo-(C_Global_Homogeneo*Velocidades_Homogeneo(:,i))-(K_Global_Homogeneo*Desplazamiento_Homogeneo(:,i)));
end 

% VISUALIZACIÓN DE LOS DESPLAZAMIENTOS, VELOCIDAD Y ACELERACIÓN DEL
% PROBLEMA HOMOGÉNEO
close all
figure
for i=1:s
subplot(3,1,1)
plot(Desplazamiento_Homogeneo(:,i))
axis([-10 length(M_Global)+10 min(min(Desplazamiento_Homogeneo)) max(max(Desplazamiento_Homogeneo))])
title('Desplazamiento')
grid on

subplot(3,1,2)
plot(Velocidades_Homogeneo(:,i))
axis([-10 length(M_Global)+10 min(min(Velocidades_Homogeneo)) max(max(Velocidades_Homogeneo))])
title('Velocidad')

subplot(3,1,3)
plot(Aceleraciones_Homogeneo(:,i))
axis([-10 length(M_Global)+10 min(min(Aceleraciones_Homogeneo)) max(max(Aceleraciones_Homogeneo))])
title('Aceleración')
drawnow
end


%% DEFINICIÓN DE LOS PARÁMETROS PARA CALCULAR LA REDUCCIÓN
Delta_M=M_Global_Homogeneo-M_Global;
Delta_C=C_Global_Homogeneo-C_Global;
Delta_K=K_Global_Homogeneo-K_Global;

% Vector de fuerzas equivalentes 
Fuerzas_equivalentes=zeros(length(Delta_M),s);
for i=1:s
    Fuerzas_equivalentes(:,i)=((Delta_M*Aceleraciones_Homogeneo(:,i))+(Delta_C*Velocidades_Homogeneo(:,i))+(Delta_K*Desplazamiento_Homogeneo(:,i)));
end 
 
% Se resuelve por medio de diferencias finitas
Delta_Ui=zeros(length(Delta_M),1);
Delta_Uim1=zeros(length(Delta_M),1);

Desplazamientos_reducido=zeros(length(Delta_M),s);
for i=1:s
    Delta_UiM1=Matriz_inversa*(((dt^2)*(Fuerzas_equivalentes(:,i)-(K_Global*Delta_Ui)))+((dt/2)*C_Global*Delta_Uim1)+(M_Global*((2*Delta_Ui)-Delta_Uim1)));
    Delta_Uim1=Delta_Ui;
    Delta_Ui=Delta_UiM1;
    Desplazamientos_reducido(:,i)=Delta_UiM1;
end 

% Se obtienen los desplazamientos del problema original (en teoría esto
% lo que se debe obtener, ya que el primer problema no se resuelve)
Desplazamiento_original=zeros(length(M_Global),s);
for i=1:s
    Desplazamiento_original(:,i)=Desplazamientos_reducido(:,i)+Desplazamiento_Homogeneo(:,i);
end

% % Se visualizan los parámetros empleados en la reducción y se comparan los problemas
% figure(1)
% subplot(1,2,1)
% spy(Delta_M)
% title('Matriz \DeltaM')
% axis([-10 length(Delta_M)+10 -10 length(Delta_M)+10])
% grid on 
% 
% subplot(1,2,2)
% spy(Delta_K)
% title('Matriz \DeltaK')
% axis([-10 length(Delta_M)+10 -10 length(Delta_M)+10])
% grid on 
% 
% sgtitle('Parámetros empleados en la reducción de dominio')
% 
% figure(2)
% subplot(2,2,1:2)
% spy(Fuerzas_equivalentes)
% title('Fuerzas equivalentes')
% axis([-10 s+10 -10 length(Delta_M)+10])
% grid on 
% 
% subplot(2,2,3:4)
% spy(Desplazamientos_reducido)
% title('Desplazamientos problema reducido')
% axis([-10 s+10 -10 length(Delta_M)+10])
% grid on 

% % Crear un objeto VideoWriter
% videoFileName = 'Solución numérica (problema homogéneo).avi'; % Nombre del archivo de video
% writerObj = VideoWriter(videoFileName, 'Motion JPEG AVI');
% % Configurar propiedades del objeto VideoWriter
% writerObj.FrameRate =100; % Velocidad de cuadros (frames por segundo)
% % Abrir el objeto VideoWriter para escribir
% open(writerObj);

% profile viewer

% close all
% figure
% for i=1:s
%     subplot(2,1,1)
%     plot(Desplazamiento(:,i),'-r')
%     axis([-20 No_nodos+20 -1 1 ])
%     % title('Problema Auxiliar (Homogéneo)')
%     title('Problema inicial')
%     grid on
% 
%     subplot(2,1,2)
%     plot(Desplazamiento_original(:,i),'-r')
%     axis([-20 No_nodos+20 -1 1 ])
%     % title('Problema Original (aproximado)')
%     title('Problema aproximado')
%     grid on
%     drawnow
%     % pause(0.1)
% % % 
% % % % %     Capturar el cuadro actual de la figura
% % % %     currentFrame = getframe(gcf);
% % % % %     Escribir el cuadro en el archivo de video
% % % %     writeVideo(writerObj, currentFrame);
% end

% % % Cerrar el objeto VideoWriter
% close(writerObj);

% Desplazamiento=zeros(length(M_Global),s);
% Fuente_desp=Fuente;
% for i=2:s-1
%     if i>2
%         Fuente_desp=Fuente*0;
%     end
%     Desplazamiento(:,i+1)=Matriz_inversa*(((dt^2)*(Fuente_desp-(K_Global*Desplazamiento(:,i))))+((dt/2)*C_Global*Desplazamiento(:,i-1))+(M_Global*((2*Desplazamiento(:,i))-Desplazamiento(:,i-1))));
% end 
% 
% Aceleraciones=zeros(length(M_Global),s);
% Fuente_ace=Fuente;
% for i=1:s
%     if i>1
%         Fuente_ace=Fuente_ace*0;
%     end 
%     Aceleraciones(:,i)=(M_Global\eye(length(M_Global)))*(Fuente_ace-(C_Global*Velocidades(:,i))-(K_Global*Desplazamiento(:,i)));
% end 

% % % Crear un objeto VideoWriter
% % videoFileName = 'Solución (problema heterogéneo).avi'; % Nombre del archivo de video
% % writerObj = VideoWriter(videoFileName, 'Motion JPEG AVI');
% % % Configurar propiedades del objeto VideoWriter
% % writerObj.FrameRate =100; % Velocidad de cuadros (frames por segundo)
% % % Abrir el objeto VideoWriter para escribir
% % open(writerObj);
% % 
% close all
% figure
% for i=1:s
%     subplot(2,1,1)
%     plot(Desplazamiento_original(:,i),'-r')
%     axis([-20 No_nodos+10 -1 1])
%     title('Problema original (aproximado)')
%     grid on
% 
%     subplot(2,1,2)
%     plot(Desplazamiento(:,i),'-r')
%     axis([-20 No_nodos+10 -1 1 ])
%     title('Problema original')
%     grid on
%     
% drawnow
% % %     Capturar el cuadro actual de la figura
% %     currentFrame = getframe(gcf);
% % %     Escribir el cuadro en el archivo de video
% %     writeVideo(writerObj, currentFrame);
% end
% % 
% % % % Cerrar el objeto VideoWriter
% % close(writerObj);


% close all

% figure
% for i=1:s
%     subplot(3,1,1)
%     plot(Desplazamiento_Homogeneo(:,i),'-r')
%     axis([-20 No_nodos+10 -1 1])
%     title('Desplazamientos (problema homogéneo)')
%     grid on
% 
%     subplot(3,1,2)
%     plot(Desplazamiento(:,i),'-r')
%     axis([-20 No_nodos+10 -1 1])
%     title('Desplazamientos (problema heterogéneo)')
%     grid on
% 
%     subplot(3,1,3)
%     plot(Fuerzas_equivalentes(:,i),'-b')
%     title('Fuerzas equivalentes')
%     axis([-20 No_nodos+10 -1 1])
%     grid on 
%     drawnow
% end

Desplazamientos=zeros(No_nodos,s);
Desplazamientos(2:end-1,1:s)=Desplazamiento_original;
% Malla para interpolar valores en la malla densa
Coordenadas_nodales=zeros(No_nodos,2);
for i=1:No_nodos
    Coordenadas_nodales(i,1)=dx*(i-1);
end 


% Ahora para organizar la información se concatenan las matrices de coordenadas
% nodales con los desplazamientos de cada uno de los nodos.
Coordenadas_desplazamientos=[Coordenadas_nodales Desplazamientos];
%Hay que recordar que las filas de la matriz corresponden a cada uno de los
%nodos del dominio, la primera columna corresponde a la coordenada nodal x,
%la segunda corresponde a la coordenada nodal en y y las tercera columna
%corresponde al desaplazmiento del nodo.


% Ya resuelto el problema heterogeneo, se genera la malla arbitraria, esta
% podrá ser de dimensión menor que la malla heterogenea (debido a esto se 
% reduce el dominio), con una menor cantidad de elementos.


%% PARÁMETROS DE LA REDUCCIÓN
Longitud_reduccion=10;
No_elementos_reduccion=1000;
dx_reduccion=Longitud_reduccion/No_elementos_reduccion;
No_nodos_reduccion=No_elementos_reduccion+1;

%Se genera la malla para el problema arbitrario
Malla_reduccion=0:dx_reduccion:Longitud_reduccion;

Coordenadas_nodales_reduccion=zeros(No_nodos_reduccion,3);
Coordenadas_nodales_reduccion(1:end,1)=Malla_reduccion;     %Este arreglo define las coordenadas de cada uno de los nodos que conforman la malla del problema heterogéneo

No_interpolaciones=length(Malla_reduccion);

%Se visualizan ambas malla
figure(1)
plot(Coordenadas_desplazamientos(:,1),Coordenadas_desplazamientos(:,2),'xr','DisplayName','Nodos del problema heterogéneo')
hold on 
plot(Coordenadas_nodales_reduccion(:,1),Coordenadas_nodales_reduccion(:,2),'xb','DisplayName','Nodos reducción')
plot(Coordenadas_desplazamientos(:,1),Coordenadas_desplazamientos(:,2),'-k','DisplayName','Dominio')
hold off
title('Superposición de mallas')
legend
grid on 
text(3,-0.7,['Número de nodos totales: ' num2str(No_nodos) ' nodos'])
text(3,-0.85,['Número de nodos (interpolación): ' num2str(No_nodos_reduccion) ' nodos'])

%Se identifica el número de elemento en donde se interpolarán cada uno de los
%valores 

% Se crea  una matriz cuya primera columna indica el número de elemento y
% las cuatro columnas siguiente contienen coordenadas, la segunda y tercera
% columna corresponden a las coordenadas del primer nodo del elemento y la
% cuarta y quinta columna son las coordenadas del segundo nodo (el nodo que cierra al elemento)
% estas coordenadas están ordenadas en forma de x y y. Las dos columnas
% restantes indican los valores de desplazamiento de cada uno de los nodos,
% valores nodales.
Valores_interpolados=zeros(No_interpolaciones,s);
for f=1:s
Matriz_coord=zeros(No_elementos,s+5);
for i=1:No_elementos
    Matriz_coord(i,1)=i;
        Matriz_coord(i,2:3)=Coordenadas_desplazamientos(i,1:2);
        Matriz_coord(i,4:5)=Coordenadas_desplazamientos(i+1,1:2);
        Matriz_coord(i,6:7)=Coordenadas_desplazamientos(i:i+1,f+2);
end

%Ahora se crea un arreglo que indique en que elementos del problema
%heterogéneo se tendrá que interpolar
Elemento_interpolacion=zeros(No_interpolaciones,1);
for i=1:No_interpolaciones  %Este ciclo identifica las coordenadas de los puntos en donde se deberá de interpolar
    Valor_a_interpolar=Malla_reduccion(1,i);                 %Coordenadas de los puntos a interpolar (Correspondiente al problema auxiliar)
    for m=1:size(Matriz_coord,1)                             %Este ciclo recorre toda la matriz del problema heterogéneo, con este ciclo se buscan las coordenadas entre las que queda el punto a interpolar
        if (Valor_a_interpolar>=Matriz_coord(m,2)) && (Valor_a_interpolar<=Matriz_coord(m,4))   %Este if identifca entre que coordenadas se encuentran las coordenadas de los putnos a interpolar
        Elemento_interpolacion(i,1)=Matriz_coord(m,1);       %Una vez identificadas las coordenadas entrre las que queda el punto a interpolar, se identifica el elemento en donde se encuetra dicho punto de interpolación 
        end 
    end
end

%Una vez identificados los elementos en donde se requiere interpolar, se
%interpolan los valores
for i=1:No_interpolaciones
    Elemento_interpolar=Elemento_interpolacion(i,1);        %Según la ubicación de los puntos de interpolación, se identifica el elemento al que pertence y de interpola sobre ese mismo elementos
    Primer_pto_interpolacion_global=Matriz_coord(Elemento_interpolar,2);
    Segundo_pto_interpolacion_global=Matriz_coord(Elemento_interpolar,4);
    %Los puntos de interpolación definidos anteriormente están en un sistema
    %global, para realizar la interpolación es necesario establecerlos, pero en
    %coordenadas locales

    Primer_pto_interpolacion_local=Primer_pto_interpolacion_global-Primer_pto_interpolacion_global;
    Segundo_pto_interpolacion_local=Segundo_pto_interpolacion_global-Primer_pto_interpolacion_global;
    Pto_interpolacion_local=Coordenadas_nodales_reduccion(i,1)-Primer_pto_interpolacion_global;
    if Pto_interpolacion_local<=0
        Pto_interpolacion_local=0;
    end 
   
    %Se realiza la intepolación lineal de los valores, con el elemento
    %definido y con los valores nodales en cada uno de ellos.
    V=Matriz_coord(Elemento_interpolar,6)*(1-(Pto_interpolacion_local/dx))+Matriz_coord(Elemento_interpolar,7)*(Pto_interpolacion_local/dx);
   
    Valores_interpolados(i,f)=V;
end 

end 
% profile viewer

%Se grafican los valores de desplazamiento para ambos problemas 
%Para el problema heterogéneo y para el problema de referencia o arbitrario
figure(2)
for i=1:s
subplot(2,1,1)
plot(Coordenadas_nodales,Desplazamientos(:,i),'-k')
% hold on
% plot(Coordenadas_nodales,Desplazamientos(:,i),'xr')
% hold off
title('Problema original')
axis([-4 Coordenadas_nodales(end,1)+4 -1 1])
grid on

subplot(2,1,2)
plot(Malla_reduccion,Valores_interpolados(:,i),'-k')
% hold on
% plot(Malla_reduccion,Valores_interpolados(:,i),'xb')
% hold off
title('Problema reducido')
axis([-4 Coordenadas_nodales(end,1)+4 -1 1])
grid on

% plot(Coordenadas_nodales,Desplazamientos(:,i),'-k')
% hold on
% plot(Malla_reduccion,Valores_interpolados(:,i),'-r')
% hold off
% axis([-1 Longitud+1 -0.75 0.75])
% grid on
drawnow
end


%Se grafican los valores de desplazamiento para ambos problemas 
%Para el problema heterogéneo y para el problema de referencia o arbitrario
figure(2)
for i=1:s
subplot(2,1,1)
plot(Coordenadas_nodales,Desplazamientos(:,i),'-k')
hold on
plot(Malla_reduccion,Valores_interpolados(:,i),'-r')
hold off
axis([-4 Coordenadas_nodales(end,1)+4 -1 1])
grid on
drawnow
end