function reconstruction = general_spectral_method(hologram, a, b, c, d, pixel_x, pixel_y)

S = size(hologram);
ny = S(1);
nx = S(2);
[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*pixel_y;
width_x = nx*pixel_x;

Chirp1 = exp(-1i*pi*(-b/a)*d*(  ((x/width_x).^2) +  ((y/width_y).^2)));
Mag = -a;
Chirp2 = exp(-1i*pi*(c/a)*d*(  ((x/width_x).^2) +  ((y/width_y).^2)));

reconstruction = (fftshift(fft2(fftshift((fftshift(fft2(fftshift(hologram)))).*Chirp1))).*Mag).*Chirp2;