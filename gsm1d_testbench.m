%% -- Example Usage of gsm1d function (initial testing phase)

clear all
fs = 2048;
t = [-(4-1/fs):1/fs:4];
%x = rectpuls(t,1);
x = sin(2*pi*5*t);
figure(1)
plot(t,x)
f = fs/length(t)*((0:length(t)-1) - length(t)/2);

alpha = 1;
a = cos(alpha*pi/2);
b = sin(alpha*pi/2);
c = -sin(alpha*pi/2);
d = cos(alpha*pi/2);
%T = transform(length(t),a,b,d,1/fs,fs/length(t));
%X = x*T;
X = gsm1d(x,1/fs, a, b, c);
%X = fft(fftshift(x));
figure(2)
plot(f, abs(X))
title('LCT Magnitude')