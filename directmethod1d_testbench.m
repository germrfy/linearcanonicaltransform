%% -- Example Usage of gsm1d function (initial testing phase)

clear all
fs = 51.2;
t = [-(5-1/fs):1/fs:5];
N = length(t)
x = rectpuls(t,1);
%x = sin(2*pi*5*t);
figure(1)
plot(t,x)

alpha = 0.5;
a = cos(alpha*pi/2);
b = sin(alpha*pi/2);
c = -sin(alpha*pi/2);
d = cos(alpha*pi/2);
y = fs*abs(b)/length(t)*((0:length(t)-1) - length(t)/2);


%T = transform(length(t),a,b,d,1/fs,fs/length(t));
%X = x*T;
X = direct_method1d(x,1/fs, a, b, c, d);
figure(2)
plot(y ,abs(X))
title('LCT Magnitude')
