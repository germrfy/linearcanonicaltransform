function lct = agsm1d(samples, T, b, c, d)
% Written by Gerard Murphy
% Function to numerically calculate the LCT of a given set of samples using
% the alternate spectral method generalisation. (Not yet Published)

nSamp = length(samples);    % get the number of samples

row = floor(-nSamp/2):floor(nSamp/2-1);   % For use to generate chirp multiplications [-N/2 -1 : N/2]

chirp1 = exp(1i*pi*((row*T).^2)*(c/d));     % first chirp

chirp2 = exp(1i*pi*((row*d/(nSamp*T)).^2)*(-b/d));     % second chirp

g = samples.*chirp1;
G = fftshift(fft(fftshift(g)));
h = G.*chirp2;
lct = fftshift(fft(fftshift(h)));
