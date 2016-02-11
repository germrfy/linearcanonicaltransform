close all;
clear all;

%change these to the paths to wherever you've stored the hologram and
%reference files.
hol_filename = 'lena_gray.tif';

period_x = 3.45e-6;
period_y = 3.45e-6;

image = imread(hol_filename);
size(image)
image = double(image(:,:,1));
size(image)

figure
imagesc(abs(image).^2)
colormap gray;
t = pi/4;

A = [1.0736 1.5290 ; -1.3651 1.5669];
B = [0.1618 -1.5123 ; -1.5477 -0.6075];
C = [-6.5937 -10.9024 ; 1.5335 -2.9495];
D = [0 10.7399 ; 1.7875 1.8244];

result = ding2DNS(image, A, B, D, period_x, period_y);

figure
imagesc(abs(result).^2)

