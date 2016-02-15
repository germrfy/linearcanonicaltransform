clear all; close all;
Lx = 6; N = 256; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2);
%f = rectpulse(x/6, N);
f = padarray(g, [0, length(g)]);
alpha = 0.3;
a = cos(alpha*pi/2);
b = sin(alpha*pi/2);
c = -sin(alpha*pi/2);
d = cos(alpha*pi/2);
%a=3;b=2;c=1;d=1;

y = abs(b)/(Tx*length(f))*((0:length(f)-1) - length(f)/2);
yg = abs(b)/(Tx*length(g))*((0:length(g)-1) - length(g)/2);

r = gsm1d(f, Tx, a, b, c);
figure; title('GSM')
subplot(411);plot(y, abs(r));
subplot(412);plot(y, angle(r));
subplot(413);plot(y, real(r));
subplot(414);plot(y, imag(r));

s = agsm1d(f, Tx, b, c, d);
figure; title('AGSM')
subplot(411);plot(y, abs(s));
subplot(412);plot(y, angle(s));
subplot(413);plot(y, real(s));
subplot(414);plot(y, imag(s));

% q = frft(f, alpha);
% q = q';
% figure; title('Ref')
% subplot(411);plot(y, abs(q));
% subplot(412);plot(y, -1*angle(q));
% subplot(413);plot(y, real(q));
% subplot(414);plot(y, imag(q));
% %% Analytic Solution (From Mathematica)
% 
% exF = exp((pi*(1i*b*d + (a*d-1)*pi)*y.^2)/(b*(b-1i*a*pi)));
% den = 2*sqrt(1i*b)*sqrt((-1i*a/b)+(1/pi))*y;
% fb = -1i*b*sqrt(1-(1i*a*pi/b))*sqrt(-(y.^2)/(b^2-1i*a*b*pi)).*(1i*erfi(pi*sqrt(-(y.^2)/(b^2-1i*a*b*pi))));
% sb = y.*(2-1i*erfi(pi*y/(b*sqrt(1-(1i*a*pi/b)))));
% 
% an = (exF.*(fb.*sb))./den;
% 
T = transform(length(g),a,b,d,1/fs,fs*(abs(b))/length(g));
q = g*T;
figure; title('Ref')
subplot(411);plot(yg, abs(q));
subplot(412);plot(yg, angle(q));
subplot(413);plot(yg, real(q));
subplot(414);plot(yg, imag(q));
% figure
% plot(y, angle(an));
