%   Written by Gerard Murphy (UCD), 31/Oct/2015

%   The following code is based on the general spectral method algorithm
%   for numerical computation of any 2D-S-DLCT as outlined in 
%   J. J. Healy and J. T. Sheridan, "Space–bandwidth ratio as a means of 
%   choosing between Fresnel and other linear canonical transform 
%   algorithms"

%	input: input 2D signal/image
%   output: output signal/image
%   period_x,period_y: sampling periods in x and y input dimensions
%   a11, b11, c11:  LCT parameters for the LCT to be performed in the x
%                   input direction
%   a22, b22, c22:  LCT paremeters for the LCT to be performed in the y
%                   input direction

function output = general_spectral_method(input, a11, a22, b11, b22, c11, c22, period_x, period_y)

S = size(input);
ny = S(1);
nx = S(2);
[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(-1i*pi*(  ((x/width_x).^2)*(b11/a11) +  ((y/width_y).^2)*(b22/a22)));
Chirp2 = exp(1i*pi*(  ((x*a11*period_x).^2)*(c11/a11) +  ((y*a22*period_y).^2)*(c22/a22)));

output = (fftshift(fft2(fftshift((fftshift(fft2(fftshift(input)))).*Chirp1)))).*Chirp2;