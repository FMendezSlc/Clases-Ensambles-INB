%%Felipe Mendez Salcido MD, MSc
% Juriquilla, Qro. Mx. 2019
function simMiniBurts = PBC_analysis(filename, channel)
% filename debe ser el archivo de Josue
% channel es el canal a analizar: 2-10
Fs = 25000; % Frecuencia de muestreo
period = 1/Fs; % periodo de muestreo en seg
%%
load(filename);
chan_1 = jCanales(channel,:); % canal de trabajo

window = round(0.5/period); % ventana de 500  ms
[up, low]= envelope(chan_1, window, 'rms'); % envolvente
up = up-median(up); % substraccion de linea basal
min_amp = median(up)+3*mad(up); % umbral de deteccion

 % Los parametros de deteccion han sido optimizados
[pks, locs] = findpeaks(up, 'MinPeakHeight', min_amp,...
    'MinPeakDistance', .5, 'MinPeakProminence', 3);
% Filtrado adicional
lpFilt = designfilt('lowpassfir', 'FilterOrder',25, ...
       'PassbandFrequency',300, 'StopbandFrequency',350, ...
'SampleRate',Fs);

Chan_filt = filtfilt(lpFilt, chan_1);

numBursts = length(pks);
allBursts = zeros(numBursts, 12500);
for ii = 1:numBursts;
    start = locs(ii)-6249;
    stop = locs(ii)+6250;
    allBursts(ii, :) = Chan_filt(start:stop);
end

corrMat = squareform(1-pdist(allBursts, 'correlation'));

%% Descomponer rafagas en mini ventanas (25ms)
mini_bursts = zeros(178*20, 625);
row_count = 1;
for ii = 1:size(allBursts, 1)
mini_start = 1;
mini_stop = 625;
for jj = 1:20
    mini_bursts(row_count, :) = allBursts(ii, (mini_start:mini_stop));
    row_count = row_count + 1;
    mini_start = mini_start + 625;
    mini_stop = mini_stop + 625;
end
end
mini_winCorr = squareform(1-pdist(mini_bursts, 'correlation'));

%diag_one = diag(ones(size(mini_winCorr)));
%mini_winCorr = mini_winCorr - diag(diag(mini_winCorr)) + diag(diag_one);

%% Crear figura; 
% trazo original con rafagas identificadas
% matriz de similitud de rafagas completas
% matris de similitud de mini ventanas
figure(1)
subplot(2,2,[1 2])
plot(chan_1);
hold on
for ii = 1: length(locs);
    start = locs(ii)-6249;
    stop = locs(ii)+6250;
    plot((start:stop), chan_1(start:stop), 'r');
end
hold off
title('Burst Detection')
subplot(2,2,3)
imagesc(corrMat)
colorbar
axis 'square'
set(gca,'YDir','normal')
title('Similarity Whole Bursts')

subplot(2, 2, 4)
imagesc(mini_winCorr);
colorbar
axis 'square'
set(gca,'YDir','normal')
title('Similarity Mini Windows')
%% Output
simMiniBurts = mini_winCorr;
