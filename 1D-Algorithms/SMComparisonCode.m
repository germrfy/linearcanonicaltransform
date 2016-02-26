clear all; close all;
Lx = 6; N = 256; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2).*rectangularPulse(x/6);
% alpha = 0.2;
% A = cos(alpha*pi/2);
% B = sin(alpha*pi/2);
% C = -sin(alpha*pi/2);
% D = cos(alpha*pi/2);
A=3;B=2;C=1;D=1;
%A=0.9;B=-1;C=0.1;D=1;

f = padarray(g, [0, round(6*length(g))]);
y = abs(B)/(Tx*length(f))*((0:length(f)-1) - length(f)/2);

%% Analytic Solution

ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

an = (ExAn.*(fb + sb))./(den);

figure;
title('Analytic Solution')
anAbs = max(abs(an));
subplot(411);plot(y, abs(an)/anAbs);title('Analytic Solution Absolute Value');
subplot(412);plot(y, angle(an));title('Analytic Solution Phase');
subplot(413);plot(y, real(an)/anAbs);title('Analytic Solution Real Part');
subplot(414);plot(y, imag(an)/anAbs);title('Analytic Solution Imaginary Part');

%% General Spectral Method

r = gsm1d(f, Tx, A, B, C);
figure; 
title('GSM')
scr = max(abs(r));
subplot(411);plot(y, abs(r)/scr);title('GSM Absolute Value');
subplot(412);plot(y, angle(r));title('GSM Phase');
subplot(413);plot(y, real(r)/scr);title('GSM Real Part');
subplot(414);plot(y, imag(r)/scr);title('GSM Imaginary Part');

%% Alternate Spectral Method

s = agsm1d(f, Tx, B, C, D);
figure; 
scs = max(abs(s));
subplot(411);plot(y, abs(s)/scs);title('AGSM Absolute Value');
subplot(412);plot(y, angle(s));title('AGSM Phase');
subplot(413);plot(y, real(s)/scs);title('AGSM Real Part');
subplot(414);plot(y, imag(s)/scs);title('AGSM Imaginary Part');

MSE_GSM = immse(r/max(abs(r)),an/max(abs(an)))
MSE_AGSM = immse(s/max(abs(s)),an/max(abs(an)))

