%% -- Example Usage of gsm1d function (initial testing phase)

clear all
fs = 128;
t = [-(1-1/fs):1/fs:1];
N = length(t)
x = rectpuls(t,1);
%x = sin(2*pi*5*t);
figure(1)
plot(t,x)

alpha = 0.1;
a = cos(alpha*pi/2);
b = sin(alpha*pi/2);
c = -sin(alpha*pi/2);
d = cos(alpha*pi/2);

% Apply Ding's Sampling theorem to samples in output domain
y = fs*abs(b)/length(t)*((0:length(t)-1) - length(t)/2);

%T = transform(length(t),a,b,d,1/fs,fs*(abs(b))/length(t));
%V = x*T;
U = direct_method1d(x,1/fs, a, b, c, d);
X = gsm1d(x,1/fs, a, b, c);
figure(2)
plot(y, abs(X))
title('LCT Magnitude SM')
% figure(3)
% plot(y, abs(V))
% title('LCT Magnitude DIR')
figure(4)
plot(y, abs(U))
title('LCT Magnitude DM')
