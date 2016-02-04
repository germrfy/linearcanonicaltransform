close all;
clear all;

%change these to the paths to wherever you've stored the hologram and
%reference files.
hol_filename = 'Hologram.tif';
ref_filename = 'ref.tif';

pixel_x = 3.45e-6;
pixel_y = 3.45e-6;

wavelength = 682.5e-9;

hol = imread(hol_filename);
ref = imread(ref_filename);

hol = double(hol(:,:,1));
ref = double(ref(:,:,1));

% Filter unwanted terms from hologram 
FFT_hol = fftshift(fft2(fftshift(hol)));
x0 = 225; y0 = 250;
R = 222;
FFT_hol = circle_filter(FFT_hol,R,x0,y0);
  %plot FT of filtered hologram
  contrastproxy = 0.1;
  MIN = min(min(abs(FFT_hol).^2)); 
  MAX = contrastproxy*max(max(abs(FFT_hol).^2));
  figure;
  imagesc(abs(FFT_hol).^2, [MIN MAX]);
  colormap gray; 
  title('Fourier Transform of Hologram');
  axis off; 
  axis equal;
hol = fftshift(ifft2(fftshift(FFT_hol)));

FFT_ref = fftshift(fft2(fftshift(ref)));
FFT_ref = circle_filter(FFT_ref,R,x0,y0);
ref = fftshift(ifft2(fftshift(FFT_ref)));

% Remove reference (hologram = f r*; reference = r; now we get f r* r = f = image)
hol = hol.*exp(i*angle(conj(ref)));

% This is the z parameter of the Fresnel transform
reconstruction_depth = -0.0025;

%recon_hol = spectral_method(hol, reconstruction_depth, pixel_x, pixel_y, wavelength);
%recon_hol = general_spectral_method(hol, 1, 1, reconstruction_depth*wavelength, reconstruction_depth*wavelength, 0, 0, pixel_x, pixel_y);
%recon_hol = general_direct_method(hol, 1, 1, reconstruction_depth*wavelength, reconstruction_depth*wavelength, 1, 1, pixel_x, pixel_y);
%recon_hol = general_algorithm3(hol, 1, 1, reconstruction_depth*wavelength, reconstruction_depth*wavelength, 0, 0, 1, 1, pixel_x, pixel_y);



contrastproxy2 = 1;
MIN = min(min(abs(recon_hol).^2)); 
MAX = contrastproxy2*max(max(abs(recon_hol).^2));
figure;
imagesc(abs(recon_hol).^2, [MIN MAX]);
colormap gray; 
title('Reconstruced Intensity');
axis off; 
axis equal;
figure;
imagesc(angle(recon_hol)); 
title('Reconstruced Phase'); 
colormap gray;
axis off; 
axis equal;