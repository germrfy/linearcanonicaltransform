%   Written by Gerard Murphy (UCD), 05/Nov/2015

%   The following code is based on the general direct method algorithm
%   for numerical computation of any 2D-S-DLCT as outlined in 
%   J. J. Healy and J. T. Sheridan, "Space–bandwidth ratio as a means of 
%   choosing between Fresnel and other linear canonical transform 
%   algorithms"

%	input: input 2D signal/image
%   output: output signal/image
%   period_x,period_y: sampling periods in x and y input dimensions
%   a11, b11, d11:  LCT parameters for the LCT to be performed in the x
%                   input direction
%   a22, b22, d22:  LCT paremeters for the LCT to be performed in the y
%                   input direction

function output = general_direct_method(input, a11, a22, b11, b22, d11, d22, period_x, period_y)

S = size(input);
ny = S(1);
nx = S(2);
[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(1i*pi*(((x*period_x).^2)*(a11/b11) +  ((y*period_y).^2)*(a22/b22)));
Chirp2 = exp(1i*pi*(((x/(width_x*b11)).^2)*(d11/b11)+(((y/(width_y*b22)).^2)*(d22/b22))));

output = fftshift(fft2(fftshift(input.*Chirp1))).*Chirp2;