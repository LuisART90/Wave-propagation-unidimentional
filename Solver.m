function [Desplazamiento,Velocidades,Aceleraciones] = Solver(dt,Fuente,Matriz_inversa,K_Global,M_Global,C_Global,s)
% Esta función resuelve la ecuación de movimiento a través del método de
% diferencias finitas de segundo orden con una aproximación central,
% también se calculan las velocidades y aceleraciones respectivas 
Desplazamiento=zeros(length(M_Global),s);
Fuente_desp=Fuente;
for i=2:s-1
    if i>2
        Fuente_desp=Fuente*0;
    end
    Desplazamiento(:,i+1)=Matriz_inversa*(((dt^2)*(Fuente_desp-(K_Global*Desplazamiento(:,i))))+((dt)*C_Global*Desplazamiento(:,i))+(M_Global*((2*Desplazamiento(:,i))-Desplazamiento(:,i-1))));
end 

% Se obtiene las velocidades
Velocidades=zeros(length(M_Global),s);
Fuente_vel=Fuente;
for i=3:s-2
    if i>3
        Fuente_vel=Fuente_vel*0;
    end
    Velocidades(:,i)=(C_Global\eye(length(M_Global)))*(((Fuente_vel-M_Global*(Desplazamiento(:,i+1)-2*Desplazamiento(:,i)+Desplazamiento(:,i-1))./(dt^2))-K_Global*Desplazamiento(:,i)));
end

Aceleraciones=zeros(length(M_Global),s);
Fuente_ace=Fuente;
for i=3:s-2
    if i>3
        Fuente_ace=Fuente_ace*0;
    end 
    Aceleraciones(:,i)=(M_Global\eye(length(M_Global)))*(Fuente_ace-(C_Global*Velocidades(:,i))-(K_Global*Desplazamiento(:,i)));
end 

end