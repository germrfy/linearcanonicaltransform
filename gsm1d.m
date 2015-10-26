function lct = gsm1d(samples, T, a, b, c)
% Written by Gerard Murphy UCD student
% Function to numerically calculate the LCT of a given set of samples using
% the spectral method.

nSamp = length(samples);    % get the number of samples
%T = (samples(nSamp) - samples(1))/nSamp;  % get spacing between samples

row = ceil(-nSamp/2):ceil(nSamp/2-1);   % For use to generate chirp multiplications

chirp1 = exp(-1i*pi*((row/T).^2)*b/a);      % first chirp
chirp2 = exp(1i*pi*((row*a*T).^2)*c/a);     % seconds chirp
                                        % final computation of LCT
lct = (fftshift(fft(fftshift(((fftshift(fft(fftshift(samples)))).*chirp1))))).*chirp2;