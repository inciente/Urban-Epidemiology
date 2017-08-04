clc
clear all

%En este programa voy a utilizar toda la información de Guadalajara para
%calcular la evolución de epidemias utilizando el modelo de contacto
%distribuidos. 

load('malla892nod.mat', 'z','mm'); 
%Triang = load('TriangFinal.csv');
%escolar = load('escolar.csv');
%adultnin = load('adultnin.csv');
 
%TriangFinal tiene X, Y, Area, Densidad y trip_attract laboral para cada
%triángulo del mapa. Escolar tiene edbásica, edmdiasuperior.

%Empezamos por obtener población total sobre y debajo de 14 años para cada
%triángulo: Usamos la dist normal de niños/total con media 0.26748 y STD
%0.082313 para ver cuántos menores viven en cada lugar
%PobFinal es X,Y, Area, PobInfantil, PobAdulta, atractINfantil,
%AtractAdulta (poblaciones absolutas)

poblacion = load('pobfinal.csv');

%Necesitamos importar las matrices de origen y destino y crear un vector de
%beta (en este caso beta será una constante /km2)

beta = poblacion(:,3)/1000000/12000;

correspond = load('correspond.csv'); 

ODSNino = load('DSNino.csv');
ODSNino = gpuArray(ODSNino);

ODSAdulto = load('DSAdulto.csv'); 
ODSAdulto = gpuArray(ODSAdulto); 
%Importamos las probabilidades de movimiento


ODINino = ODSNino*0.3; centromedico = [115, 189];
%ODINino(:,centromedico) = ODINino(:,centromedico) + 0.03;
ODINino = gpuArray(ODINino); 

ODIAdulto = ODSAdulto*0.3; 
%ODIAdulto(:,centromedico) = ODIAdulto(:,centromedico) + 0.03;
ODIAdulto = gpuArray(ODIAdulto); 

%Vamos a definir adónde van al hospital las personas

for k = 1:1580
   ODINino(k,correspond(k)) = ODINino(k,correspond(k)) + 0.06;
   ODIAdulto(k, correspond(k)) = ODIAdulto(k,correspond(k)) + 0.06;
end


tstep = 75;

resultadoN = NaN(1580); 
resultadoA = NaN(1580); 

%parpool(2)

tic

parfor k = 1:1580
    
    Inin = gpuArray(zeros(1580,1)) ; Iadu = gpuArray(zeros(1580,1));
    Iadu(k) = 5; %Condicion inicial
    %Iadun = Iadu;
    
    k
    for t=1:tstep  
        viaI = ODIAdulto'*Iadu + ODINino'*Inin;
        %viaInin = ODINino'*Inin;
        %viaIadu = ODIAdulto'*Iadu;
        
        %Calculamos la poblacion sana en cada parte de la ciudad
        ninsanos = poblacion(:,4) - Inin(:); 
        adusanos = poblacion(:, 5) - Iadu(:);
        
        %Ahora calculamos cada renglon de ODS por S para tener valor esperado
        %de sanos que van de x a y'.
        viaSadu = bsxfun(@times, ODSAdulto, adusanos);
        viaSnin = bsxfun(@times, ODSNino, ninsanos); 
        
        %Aqui ya vemos cuando interactucan con enfermos. Empezamos por
        %distinguir transmisión de nino a adulto como nin_adu, de adulto a
        %nino como adu_nin. 
        %nin_adu = viaSadu*(beta.*viaInin); 
        %nin_nin = viaSnin*(beta.*viaInin);
        
        %adu_nin = viaSnin*(beta.*viaIadu); 
        %adu_adu = viaSadu*(beta.*viaIadu); 
        
        Iadu = Iadu + (viaSadu*(beta.*viaI));
        Inin = Inin + (viaSnin*(beta.*viaI)); 
        
    end
    
    resultadoA(:,k) = gather(Iadu); 
    resultadoN(:,k) = gather(Inin); 
end


toc