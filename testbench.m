%% -- Example Usage of gsm1d function (initial testing phase)

clear all
fs = 512;
t = [-(0.1-1/fs):1/fs:0.1];
%t = [-(2-1/fs):1/fs:2];
N = length(t)
%x = rectpuls(t,1);
x = exp(-(t/0.02).^2).*cos(100*pi*t).*exp(3*1i*pi*t); % Function Taken from Liang's work
figure(1)
plot(t,abs(x))

alpha = 0.125;
a = cos(alpha*pi/2);
b = sin(alpha*pi/2);
c = -sin(alpha*pi/2);
d = cos(alpha*pi/2);

% Apply Ding's Sampling theorem to samples in output domain
y = fs*abs(b)/length(t)*((0:length(t)-1) - length(t)/2);

% T = transform(length(t),a,b,d,1/fs,fs*(abs(b))/length(t));
% V = x*T;
U = direct_method1d(x,1/fs, a, b, c, d);
X = algorithm3(x,1/fs, a, b, c, d);
figure(2)
plot(y, abs(X))
title('LCT Magnitude SM')
% figure(3)
% plot(y, abs(V))
% title('LCT Magnitude N-Squared')
figure(4)
plot(y, abs(U))
title('LCT Magnitude DM')
% figure(5)
% plot(y, (abs(U))./(abs(V)))
% title('LCT Magnitude DM / N-Squared')
