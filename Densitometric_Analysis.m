clear all, clc, close all

img = imread('E:\Timm\300419b_Male_KO_5der.tif');

scale = round(0.625^2, 39);
[nrows, ncols, ~] = size(img);
imshow(img);
user_rotate = 'Rotate?';


if input(user_rotate) == 1
    
    grados = input('Grados:');
    rotat =imrotate(img, grados);
    imshow(rotat); % Imagen original
    
    
    h = imrect;
    position = wait(h); % Seleccionar ROI
    croppedImage = imcrop(rotat, position);
else
    imshow(img); % Imagen original
    h = imrect;
    position = wait(h); % Seleccionar ROI
    croppedImage = imcrop(img, position);
end
figure;
imshow(croppedImage);
bwimg = rgb2gray(croppedImage); % Convertir a grayscale
imshow(bwimg)
% fondo = imrect; % lineas para restar fondo como en el western
% pos2 = wait(fondo);
% fondoSel = imcrop(bwimg, pos2);
% rev_img = 255-bwimg -(255-mean(fondoSel(:))); % Restar intensidad de la imagen
fondo = input('Fondo: ');
rev_img = 255-bwimg - fondo; % Restar intensidad de la imagen
binaryImage =(rev_img > 100);
imshow(binaryImage);
colormap(gray)
Medidas = regionprops(binaryImage);
shg

AreaTimm = sum([Medidas.Area]); % Suma de las areas que revasan umbral.
AreaMicrons = AreaTimm*scale


% figure(4)
% rectangle('Position', [0,0, nrows, ncols], 'Facecolor', 'w'), hold on
% rectangle('Position', [0,0, sqrt(AreaTimm), sqrt(AreaTimm)], 'Facecolor', 'k')
% axis('tight')
% title('Area de la region seleccionada vs area total')
