function output = ding2DNS(input, A, B, D, period_x, period_y)

C1 = inv(B)*A;      % Ignore warnings as only 2x2 matrices
C2 = D*inv(B);


S = size(input);
ny = S(1);
nx = S(2);

[x, y]= meshgrid(ceil(-nx/2):ceil(nx/2-1),ceil(-ny/2):ceil(ny/2-1));

width_y = ny*period_y;
width_x = nx*period_x;

Chirp1 = exp(-1i*pi*(  ((x*period_x).^2)*(C1(1,1)) + (x*y*period_x*period_y)*(C1(1,2) + C1(2,1)) + ((y*period_y).^2)*(C1(2,2))) );

%   Affine Transformation
Tm = [[A [0,0]'] ; [0,0,1]];    % Affine Transfermation Matrix
%Imwarp(input, Tm);

Chirp2 = exp(1i*pi*(  ((x/width_x).^2)*(C2(1,1)) + ((x*y/(width_x*width_y))*(C2(1,2) + C2(2,1))) + ((y/width_y).^2)*(C2(2,2)))); % Check this


output = Imwarp(fftshift(fft2(fftshift(input.*Chirp1))), Tm).*Chirp2;