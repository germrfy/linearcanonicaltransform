function lct = direct_method1d(samples, T, a, b, c, d)
% Written by Gerard Murphy - UCD ME student
% Function to numerically calculate the LCT of a given set of samples using
% the spectral method.

nSamp = length(samples);    % get the number of samples

row = ceil(-nSamp/2):ceil(nSamp/2-1);   % For use to generate chirp multiplications

chirp1 = exp(1i*pi*((row*T).^2)*a/b);     % first chirp
chirp2 = exp(1i*pi*((row*(1/(nSamp*T*b))).^2)*d/b);     % second chirp
                                        % final computation of LCT
lct = fftshift(fft(fftshift(samples.*chirp1))).*chirp2;