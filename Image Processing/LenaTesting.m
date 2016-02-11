close all;
clear all;

%change these to the paths to wherever you've stored the hologram and
%reference files.
hol_filename = 'lena_gray.tif';

image = imread(hol_filename);
image = double(image(:,:,1));

figure
imagesc(abs(image).^2)
colormap gray;
t = pi/4;
%tform = affine2d([cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1]);
t = affine2d([2 0 0;0 0.5 0 ; 0 0 1]);
Tm = [2 0 0; 0 0.5 0; 0 0 1];    % Affine Transfermation Matrix

tform = maketform('affine', Tm);


[J,cdata,rdata] = imtransform(image,tform);
J = imwarp(image,t, 'linear');

figure
imagesc(abs(J).^2)
colormap gray;

%image = double(image1(:,:,1));

%  I = image;
%  class_of_I = class(I);
%  [x, y] = meshgrid(1:512);
%  [xi, yi] = meshgrid(1:0.1:512);
%  class_of_I
%  New_Image = cast(interp2(x,y,double(I),xi,yi,'linear'),class_of_I);

% x = fftshift(fft2(fftshift(image)));
% MIN = min(min(abs(x).^2)); 
% MAX = 0.1*max(max(abs(x).^2));
% figure;
% imagesc(abs(x).^2, [MIN MAX]);
% title('Image');
% axis off; 
% axis equal;
% 
% x = fftshift(ifft2(fftshift(image)));
% figure;
% imagesc(real(x), [MIN MAX]);
% title('Image');
% axis off; 
% axis equal;