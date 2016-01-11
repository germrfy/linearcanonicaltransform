function reconstruction = general_spectral_method(hologram, a11, a22, b11, b22, c11, c22, period_x, period_y)

S = size(hologram);
ny = S(1);
nx = S(2);
[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(-1i*pi*(  ((x/width_x).^2)*(b11/a11) +  ((y/width_y).^2)*(b22/a22)));
Chirp2 = exp(1i*pi*(  ((x*a11*period_x).^2)*(c11/a11) +  ((y*a22*period_y).^2)*(c22/a22)));

reconstruction = (fftshift(fft2(fftshift((fftshift(fft2(fftshift(hologram)))).*Chirp1)))).*Chirp2;