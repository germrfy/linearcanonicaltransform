function reconstruction = general_direct_method(hologram, a11, a22, b11, b22, d11, d22, period_x, period_y)

S = size(hologram);
ny = S(1);
nx = S(2);
[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(1i*pi*(((x*period_x).^2)*(a11/b11) +  ((y*period_y).^2)*(a22/b22)));
Chirp2 = exp(1i*pi*(((x/period_x*b11).^2)*(d11/b11)+(((y/period_y*b22).^2)*(d22/b22))));

reconstruction = fftshift(fft2(fftshift(hologram.*Chirp1))).*Chirp2;