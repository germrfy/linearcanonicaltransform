function lct = gsm1d(samples, T, a, b, c)
% Written by Gerard Murphy
% Function to numerically calculate the LCT of a given set of samples using
% the spectral method.
%samples = samples';
nSamp = length(samples);    % get the number of samples

row = floor(-nSamp/2):floor(nSamp/2-1);  % For use to generate chirp multiplications [-N/2 -1 : N/2]

chirp1 = exp(1i*pi*((row/(nSamp*T)).^2)*(-b/a));     % first chirp

chirp2 = exp(1i*pi*((row*T*a).^2)*(c/a));     % second chirp


F = fftshift(fft(fftshift(samples)));
g = F.*chirp1;
G = fftshift(fft(fftshift(g)));
lct = G.*chirp2;

