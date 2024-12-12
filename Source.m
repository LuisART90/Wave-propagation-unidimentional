function [Fuente] = Source(Distancia_fuente,dx,Frecuencia_simulacion,Longitud,dt,Nodos_fijos)
% Esta función construye la fuente empleada para la simulación numérica de
% la ecuación de onda, debido a que se considera que la fuente es impulsiva
% la fuente se construye como una función Gaussiana
% Función Gaussiana 
mu=Distancia_fuente;
sigma=dx;
if Frecuencia_simulacion>1
x=0:dx:Longitud;
elseif Frecuencia_simulacion<=1
    x=0:dx:Longitud;
end
Gaussian=((((sigma*sqrt(2*pi)))*exp(-(((x-mu).^2)/(2*sigma*sigma))))/(dt))';
Fuente=Gaussian;
Fuente(Nodos_fijos,:)=[];
end