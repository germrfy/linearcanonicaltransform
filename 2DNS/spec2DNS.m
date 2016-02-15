function [output, samx, samy] = spec2DNS(input, A, B, C, period_x, period_y)

C1 = -inv(A)*B;      % Ignore warnings as only 2x2 matrices
C2 = C*inv(A);


S = size(input);
ny = S(1);
nx = S(2);

[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(-1i*pi*(  ((x/width_x).^2)*(C1(1,1)) + ((x.*y)/(width_x*width_y))*(C1(1,2)+ C1(2,1)) + ((y/width_y).^2)*(C1(2,2))) );

%   Affine Transformation
Tm = affine2d([[A [0,0]'] ; [0,0,1]]);    % Affine Transfermation Matrix

FT1 = fftshift(fft2(fftshift(input))).*Chirp1;
J = imwarp(fftshift(fft2(fftshift(FT1))), Tm);

S = size(J);
ny = S(1);
nx = S(2);

d = inv((-A)');

newSampx = d(1,1)*(period_x) + d(2,1)*(period_y);
newSampy = d(1,2)*(period_x) + d(2,2)*(period_y);

[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

Chirp2 = exp(1i*pi*(  ((x*(newSampx)).^2)*(C2(1,1)) + (((x.*y)*(newSampy*newSampx))*(C2(1,2)+ C2(2,1))) + ((y*(newSampy)).^2)*(C2(2,2)))); % Check this

output = J.*Chirp2;
samx = newSampx;
samy = newSampx;

%output = affine2d(fftshift(fft2(fftshift(input.*Chirp1))), Tm).*Chirp2;
%output = imwarp(fftshift(fft2(fftshift(input.*Chirp1))), Tm).*Chirp2;
