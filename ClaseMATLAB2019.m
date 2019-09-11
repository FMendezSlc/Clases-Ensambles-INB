clear all, clc, close all

Re = -2:.005:2;
Im = -2:.005:2;
numCiclos = 1000;
ren = 0;

for R = Re
    ren = ren +1;
    col = 0;
    for I = Im
        col = col +1;
        contador = 1;
        c = R+I*i;
        z = c;
        while (~isinf(abs(real(z))) || ~isinf(abs(imag(z))))...
                && contador < numCiclos
            contador = contador +1;
            z = z^2+c;
        end
        if contador ==numCiclos
            pertenece(ren, col) = 1;
        else
            pertenece(ren, col) = 0;
        end
    end
end

imagesc(pertenece')

return
%Logistic growth. Ecuacion recursiva

N = 11;
r  = 1; 
k = 200;

for r= 2.8:.01:4.5
for t = 1:5000
N(t+1) = N(t)*r * (1-N(t)/k);
end

unicos = unique(N(end-150:end));
h = plot(r, unicos, 'k.-');, hold on
set(h, 'markersize', .5);
drawnow
pause(0.05)
end


return
plot(N, '*-')
xlabel('tiempo'), ylabel('# Individuos')



return
Fs = 44100; % CD quality
t = 0:1/Fs:1;
y = chirp(t, 200, 1, 2000, 'quadratic', [], 'convex');
%sound(y, Fs)
spectrogram(y,1000, 800,200:10:2000, Fs,  'yaxis')

return 
% Utilizando Spike sorting, usar PCA para reduccion de dimensiones
% Graficar los primeros 3 PCAs
% Usar K means para encontrar clusters
load('SpikeSorting')

espigas  = session.wf;
plot(espigas(30:50,:)')

return 

datosIris = iris_dataset;
datosIris = datosIris';

figure(1)
subplot(311)
stem(datosIris(:,1)), xlabel('# Flor'), ylabel('Sepal length')
subplot(312)
plot(datosIris(:,1), datosIris(:,2), 'o')
xlabel('Sepal length'), ylabel('Sepal width')
subplot(313)
plot3(datosIris(:,1), datosIris(:,2), datosIris(:,3), 'o')
xlabel('Sepal length'), ylabel('Sepal width'), zlabel('Petal length')

[coeff, datosPCA, latent] = pca(datosIris);

figure(2), hold on
plot(datosPCA(:,1), datosPCA(:,2), 'o')
xlabel('PCA1'), ylabel('PCA2'), zlabel('PCA3')

idx = kmeans(datosPCA(:,1:3), 3);
colores = [138 43 226; 64 224 208; 205 92 92]./255;
for n = 1:size(datosPCA, 1)
    plot(datosPCA(n, 1), datosPCA(n,2), 'o', 'color', colores(idx(n),:))
end
shg
return
% PCAs
n = 1000;   % numero de observaciones
rho = 0.6;  % Coeficiente de correlacion
Q = ones(2) * rho;
Q(1,1) = 1;
Q(2,2) = 1;
Q = sqrtm(Q);

rnd = randn(2, n);
rndC = Q*rnd;
datos = rndC';

subplot(211)
plot(datos(:,1), datos(:,2), 'o'), title('Datos crudos'), xlabel('X')
ylabel('Y'), hold on

Cov = cov(datos)  % Matriz de covarianza
Var = var(datos)
miVar = sum((datos(:,1)-mean(datos(:,1))).^2)./(n-1)
miCov = sum((datos(:,1)-mean(datos(:,1))).*...
            (datos(:,2)-mean(datos(:,2))))./(n-1)

Sigma = cov(datos)
[V D] = eig(Sigma)

plot([0 V(1,1)], [0 V(1,2)], 'g', 'Linewidth', 3)
plot([0 V(2,1)], [0 V(2,2)], 'r', 'Linewidth', 3)

coeff = fliplr(V);
datosPCA = datos * coeff;

subplot(212)
plot(datosPCA(:,1), datosPCA(:,2), 'o')
title('Datos reorientados'), xlabel('PCA1'), ylabel('PCA2')

shg

return
% Tarea: Crear 2 senoidales Frecs 62 y 77
% Agregarle ruido 1/f (colored noise, pink noise)
% Calcular el espectro de frecuencia de la senial resultante
% A random walk

azar = randn(1000, 1);
randWalk = cumsum(azar);
plot(azar), hold on
plot(randWalk)


return 

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t) + .1*sin(2*pi*215*t);

X = S + 2*randn(size(t));
fc = 100;
% Genera los coeficientes


subplot(211)
plot(1000*t(1:100),X(1:100))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(212)
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


return
% Fourier Transform
senial = [0 -1 -2 -3 -4 -5 -6 -7 -8 0];
x = 0:.05:2*pi;
fftSenial = fft(senial);
N = length(fftSenial);
ejeTiempo = linspace(0,2*pi,N+1);

% Magnitud

magnitud = abs(fftSenial);
angulo = angle(fftSenial);

for s = 1:N  % Crear y sumar senoidales
    amplitud = abs(fftSenial(s))/N;
    frec = s-1;
    fase = angle(fftSenial(s));
    y = amplitud*cos(x*frec+fase);
    ySenial = amplitud*cos(ejeTiempo*frec+fase);
    cosenos(s,:) = y;
    cosenosSenial(s,:) = ySenial;
end
figure(1), subplot(211)
plot(cosenos')
sumaCos = sum(cosenos, 1); % suma por columnas
sumaCosSen = sum(cosenosSenial, 1);
subplot(212)
plot(x, sumaCos, 'b-'), hold on
plot(ejeTiempo(1:end-1), sumaCosSen(1:end-1), 'r*')


return
% Example sum of senoidals

x = 0:.05:2*pi;
frec = 3;
fase = pi/3; 
y1 = 1*cos(x*frec+fase);
y2 = 2*cos(x*2.3*frec+fase);

figure(1), subplot(211)
plot(x, y1, 'r-', x, y2, 'k-')
subplot(212)
plot(x, y1+y2, 'b-')

return 
%Tarea tiempos de reaccion en la maniana y tarde
% Calcular la diferencia entre los promedios 
% Estimar la probabilidad de que dicha diferecnai se deba al asar


datosTiemposReaccion2019
whos
centrosDeBin = 0.05:.05:1.2;
% Data Conditioning
selectEnsayosH = RThombres<1;
selectEnsayosM = RTmujeres<1;
RThombres = RThombres(selectEnsayosH)
RTmujeres = RTmujeres(selectEnsayosM)
selectEnsayosH = RThombres>.1;
selectEnsayosM = RTmujeres>.1;
RThombres = RThombres(selectEnsayosH)
RTmujeres = RTmujeres(selectEnsayosM)

figure(1)
subplot(211)
hist(RThombres, centrosDeBin), title('RTh')
subplot(212)
hist(RTmujeres, centrosDeBin), title('RTm')

diferenciaRT = abs(mean(RThombres)-mean(RTmujeres))

nDatosH = length(RThombres);
nDatosM = length(RTmujeres);
nRemuestras = 2000
for r = 1:nRemuestras
    remuestraH = RThombres(randi(nDatosH, nDatosH, 1));
    remuestraM = RTmujeres(randi(nDatosM, nDatosM, 1));
    promedioRH = mean(remuestraH);
    promedioRM = mean(remuestraM);
    differenciaPromedios(r) = promedioRH-promedioRM;
end
%[mean(remuestraH) mean(remuestraM)]
hipotesisNula = differenciaPromedios-mean(differenciaPromedios);
p = sum(diferenciaRT<hipotesisNula)/nRemuestras;
if p == 0;
    p = 1/nRemuestras;
end
figure(2)
hist(hipotesisNula, -.3:.01:.3), hold on
disp(['Significancia de la diferencia es ' num2str(p)])
plot([diferenciaRT diferenciaRT], [0 350], 'r-')
return
numEnsayos = 20;

RT = nan(1, numEnsayos);
for n = 1:numEnsayos
    preStim = rand*2+1;
    pause(preStim)
    disp('Responde')
    tic
    pause
    RT(n) = toc;
    disp(['Ensayo No: ' num2str(n)])
    disp(['RT: ' num2str(RT(n))])
    disp('Listo?')
    pause
end
save('tiemposReaccionFMS', 'RT')

return
% tomar imagen, mostrar imagen original
% 3 subplots de intensidad por colores por canal RGB

%hipp = imread('hipp_cajal.jpg')
%image(hipp)
%axis equal
%shg

filtro  = [.11]';
filtro = filtro/sum(filtro);
hippFilt = convn(hipp, filtro, 'same');
figure(1)
subplot(311)
image(hipp)
subplot(312)
image(hippFilt)
subplot(313)
promedio = mean(hipp(:));
hipp2 = double(hipp);
hippCont = (hipp2-promedio)*2+promedio;
image(hippCont)

return
clear all, clc
a = [1 2 3; 4 5 6;  7 8 9];
a = 100-spiral(10);
mapaC = gray(length(a(:)));
surf(a)
colormap(mapaC)
colorbar
shg

return
a = [1 2 3; 4 5 6;  7 8 9];

mapColores = [...
    1 0 0
    0 1 0
    .5 .5 .5
    1 1 1 
    0 0 0
    0 1 1
    1 1 0
    .5 .9 .5
    .5 .3 .9];
image(a)
%surf(a);
colormap(mapColores)
shg

return
senial = zeros(1, 100);
senial(45:55) = [1 0 0 0 0 0 0 0 0 0 1];
senial(70:100) = ones(1, 31);
filtro = [16 8 4 2 1];
filtro = filtro/sum(filtro);


figure(1)
senialFiltrada = conv(senial, filtro);
stem(senial, 'k.-'), hold on
stem(senialFiltrada, 'r.-', 'linewidth', 2), hold on
shg

return
hipp = imread('hipp_cajal.jpg')
image(hipp)
axis equal
shg




return
datos = [...
    0 .01
    1 .09
    2 .21
    3 .29
    4 .35
    5 .39
    6 .40
    7 .42
    8 .41
    9 .415 ];
S = datos(:,1);
velReac = datos(:,2)';
figure(1)
h = plot(S, velReac, '*');
hold on
set(h, 'markersize', 5, 'markerfacecolor', 'k')
xlabel('[S]'), ylabel('Vel Reaccion')


Vmax = .45;
Km = 1.7;
Sajuste = 0:.1:.9;
velReacAjust = michaelisMenten(Sajuste, Vmax, Km);
plot(Sajuste, velReacAjust, 'r-')
shg

errorDeAjuste = sum((velReacAjust-velReac).^2);
contRen  = 0;

for Vmax = 0.01:.05:1
    contRen = contRen +1;
    contCol = 0;
    for Km = .01:.05:2
        contCol = contCol+1;
        velReacAjust = michaelisMenten(Sajuste, Vmax, Km);
        errorDeAjuste = sum((velReacAjust-velReac).^2);
        matrizDeError(contRen, contCol) = errorDeAjuste;
    end
end
figure(2)
surf(.01:.05:2, 0.01:.05:1, matrizDeError, 'EdgeColor', 'none');
xlabel('Km'), ylabel('Vmax'), zlabel('Error');

errorMinimo = min(matrizDeError(:))
[ren col] = find(matrizDeError==errorMinimo)
%%
return
% Michaelis-Menten equation

S = 0:.1:30;
Vmax = 3;
Km = 2;
velReact = (Vmax.*S)./(Km+S);
plot(S, velReact, '+--', 'linewidth', .5)
xlabel('[S]'), ylabel('Vel Reaccion')
shg
return
% Generar la funcion miCirculo(radio, centro, x)
% Graficar un circulo con radio 3, centro (a, b)


% x.^2+y.^2=1
radio = 1;
x = 0:pi/100:pi;
x = cos(x);
y1 = sqrt(radio-x.^2); %funcion del circulo unitario
y2 = -y1;
y = [y1 y2];
x = [x fliplr(x)];
plot(x, y, '*-')


return
%
m = 2
b = 1;
x = -10:10
y = miLinea(x, m, b);
plot(x, y, '*-')

return
x = 0:.1:5;
y = exp(x);
plot(x,y, '*-')

return
format bank
dinero = 1;
plazos = 365*24;
for numReinver = 1:plazos
    ahorro(numReinver)= dinero*(1 + 1/numReinver)^numReinver;
end
numReinver = 1:plazos;
plot(numReinver, ahorro, '*-')
xlabel('Numero de reinversiones')
ylabel('Capital')


return
primos = [];
x = 1;
tiempos = [];
nPrimos = [];
for numPrimos = 50:500:2e3
    tic
    while length(primos)< numPrimos
        x = x +1;
        esPrimo = 1;
        for divisor = 2:x-1
            residuo = rem(x,divisor);
            if residuo == 0
                esPrimo = 0;
            end
        end
        if esPrimo
            primos(length(primos)+1) = x;
        end
    end
    tiempos(length(tiempos)+1) = toc;
    nPrimos(length(nPrimos)+1) = numPrimos;
end
plot(nPrimos, tiempos, '*-')

return 
Fs = 4e3;
t = 0:1/Fs:2;
y = chirp(t, 0, 250);
sound(y, Fs)
plot(t, y)



return
% tarea usar un for loop para encontrar los primeros 50 numeros primos
% usar rem() en un for loop con if
% time it (50)
% time ti 100, 200, ... 1000 vector

fib = [1 1];
for x = 3:100
   fib(x) = fib(x-1)+fib(x-2)
   fib(x)/fib(x-1)
end

return
muestra = randn;
if muestra > 0
    display('yep')
else
    display('nop')
end


return
y = randn(1,1000)*10;
hist(y, 100), shg
desEst = std(y)
prom = mean(y)

return
Fs = 44000*3;
y = randn(1,Fs);
sound(y, Fs)


return
numeros = randn(1,1e3);
hist(numeros), shg

return
numeros = rand(1,1e5)*10;
histogram(numeros,0:.1:10)

return
Fs = 44000;
t = 0:1/Fs:1;
Frec = 2e3;
y = sin(t*2*pi*Frec);
y2 = sin(t*2*pi*Frec*2);
%plot(y)
%shg
sound(y+y2, Fs)

return
% Primer script de MATLAB
clc
disp('Hello World')

t = 0:.01:1; % tiempo
freq = 8; % Hz
amp = 9; % a.u.
phase = 2*pi;
y = sin(2*pi*t*freq + phase)*amp;
h = plot(t, y, '*b-');
xlabel('tiempo (s)');
ylabel('seno [a.u.]')
set(gca, 'fontsize', 15)
shg