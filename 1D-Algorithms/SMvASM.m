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
figure; 
title('GSM')
scr = max(abs(r));
subplot(411);plot(y, abs(r)/scr);title('GSM Absolute Value');
subplot(412);plot(y, angle(r));title('GSM Phase');
subplot(413);plot(y, real(r)/scr);title('GSM Real Part');
subplot(414);plot(y, imag(r)/scr);title('GSM Imaginary Part');

s = agsm1d(f, Tx, b, c, d);
figure; 
scs = max(abs(s));
subplot(411);plot(y, abs(s)/scs);title('AGSM Absolute Value');
subplot(412);plot(y, angle(s));title('AGSM Phase');
subplot(413);plot(y, real(s)/scs);title('AGSM Real Part');
subplot(414);plot(y, imag(s)/scs);title('AGSM Imaginary Part');

T = transform(length(g),a,b,d,1/fs,fs*(abs(b))/length(g));
q = g*T;
qI = interp1(yg,q,y,'spline');
figure; 
scq = max(abs(qI));
subplot(411);plot(y, abs(qI)/scq);title('Ref Absolute Value');
subplot(412);plot(y, angle(qI));title('Ref Phase');
subplot(413);plot(y, real(qI)/scq);title('Ref Real Part');
subplot(414);plot(y, imag(qI)/scq);title('Ref Imaginary Part');


MSE_GSM = immse(abs(r)/max(abs(r)),abs(qI)/max(abs(qI)))
MSE_AGSM = immse(abs(s)/max(abs(s)),abs(qI)/max(abs(qI)))
