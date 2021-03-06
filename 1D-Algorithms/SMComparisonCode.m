clear all; close all;
Lx = 6;Ly = 6; N = 32; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2).*rectangularPulse(x/6);
figure
plot(x, g);
title('Gaussian Input')
xlabel('x');
ylabel('f(x)')
% alpha = 0.2;
% A = cos(alpha*pi/2);
% B = sin(alpha*pi/2);
% C = -sin(alpha*pi/2);
% D = cos(alpha*pi/2);
% A=3;B=2;C=1;D=1;
% A=0.9;B=1;C=-0.1;D=1;     %%PAPER
 A=0.5;B=4;C=-0.125;D=1;   %%PAPER
 %A=1;B=2;C=2.5;D=6;        %%PAPER
 %A=0.1;B=0.9;C=-1;D=1;
f = padarray(g, [0, 4*round(length(g))]);
yGSM = abs(A)*Tx*((0:length(f)-1) - length(f)/2);       %output vector GSM
yAGSM = (Tx/abs(D))*((0:length(f)-1) - length(f)/2);    %output vector AGSM
yF= ((2*pi*fs)/(length(f)))*((0:length(f)-1) - length(f)/2); %output for FFT
yDM= (abs(B)/(length(f)*Tx))*((0:length(f)-1) - length(f)/2);

y = [yGSM; yAGSM; yDM];

%% Analytic Solution [2 analytic solns for GSM, AGSM (sample spacing)]

ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

an = (ExAn.*(fb + sb))./(den);
% plot one of the analytic solutions
figure;
anAbs = max(abs(an(1,:)));
subplot(411);plot(y(1,:), abs(an(1,:))/anAbs);title('Analytic Solution Absolute Value');
subplot(412);plot(y(1,:), angle(an(1,:)));title('Analytic Solution Phase');
subplot(413);plot(y(1,:), real(an(1,:))/anAbs);title('Analytic Solution Real Part');
subplot(414);plot(y(1,:), imag(an(1,:))/anAbs);title('Analytic Solution Imaginary Part');
figure;
subplot(211);plot(y(1,:), abs(an(1,:)));title('Analytic Solution Magnitude');xlabel('y');ylabel('|F(y)|', 'Interpreter', 'tex')
subplot(212);plot(y(1,:), angle(an(1,:)));title('Analytic Solution Phase');xlabel('y');ylabel('\angle F(y)', 'Interpreter', 'tex')

%% FT Analytic Solution

FT = (exp(-(yF.^2/4)).*(erfz(3-(1i*yF/2))+erfz(3+(1i*yF/2))))/(2*sqrt(2));
figure;
subplot(211);plot(yF, abs(FT));title('FT Analytic Solution Magnitude');xlabel('y');ylabel('|F(y)|', 'Interpreter', 'tex')
subplot(212);plot(yF, angle(FT));title('FT Analytic Solution Phase');xlabel('y');ylabel('\angle F(y)', 'Interpreter', 'tex')

%% Direct Method

dm = direct_method1d(f, Tx, A,B,D);
figure; 
scdm = max(abs(dm));
subplot(411);plot(y(3,:), abs(dm)/scdm, 'ob');hold on;plot(y(3,:), abs(an(3,:))/anAbs, '.r');title('DM Absolute Value Vs Analytic Solution');
subplot(412);plot(y(3,:), angle(dm), 'ob');hold on;plot(y(3,:), angle(an(3,:)), '.r');title('DM Phase Vs Analytic Solution');
subplot(413);plot(y(3,:), real(dm)/scdm, 'ob');hold on;plot(y(3,:), real(an(3,:))/anAbs, '.r');title('DM Real Part Vs Analytic Solution');
subplot(414);plot(y(3,:), imag(dm)/scdm, 'ob');hold on; plot(y(3,:), imag(an(3,:))/anAbs, '.r');title('DM Imaginary Part Vs Analytic Solution');


%% General Spectral Method

r = gsm1d(f, Tx, A, B, C);
figure; 
scr = max(abs(r));
subplot(411);plot(y(1,:), abs(r)/scr, 'ob');hold on;plot(y(1,:), abs(an(1,:))/anAbs, '.r');title('GSM Absolute Value Vs Analytic Solution');
subplot(412);plot(y(1,:), angle(r), 'ob');hold on;plot(y(1,:), angle(an(1,:)), '.r');title('GSM Phase Vs Analytic Solution');
subplot(413);plot(y(1,:), real(r)/scr, 'ob');hold on;plot(y(1,:), real(an(1,:))/anAbs, '.r');title('GSM Real Part Vs Analytic Solution');
subplot(414);plot(y(1,:), imag(r)/scr, 'ob');hold on; plot(y(1,:), imag(an(1,:))/anAbs, '.r');title('GSM Imaginary Part Vs Analytic Solution');

%% Alternate Spectral Method

s = agsm1d(f, Tx, B, C, D);
figure;
anAbs = max(abs(an(2,:)));
scs = max(abs(s));
subplot(411);plot(y(2,:), abs(s)/scs, 'ob');hold on;plot(y(2,:), abs(an(2,:))/anAbs, '.r');title('AGSM Absolute Value Vs Analytic Solution');
subplot(412);plot(y(2,:), angle(s), 'ob');hold on;plot(y(2,:), angle(an(2,:)), '.r');title('AGSM Phase Vs Analytic Solution');
subplot(413);plot(y(2,:), real(s)/scs, 'ob');hold on;plot(y(2,:), real(an(2,:))/anAbs, '.r');title('AGSM Real Part Vs Analytic Solution');
subplot(414);plot(y(2,:), imag(s)/scs, 'ob');hold on; plot(y(2,:), imag(an(2,:))/anAbs, '.r');title('AGSM Imaginary Part Vs Analytic Solution');

FF = fftshift(fft(fftshift(f)));
figure;
anAbs = max(abs(FT));
scF = max(abs(FF));
subplot(411);plot(yF, abs(FF)/scF, 'ob');hold on;plot(yF, abs(FT)/anAbs, '.r');title('FFT Absolute Value Vs Analytic Solution');
subplot(412);plot(yF, angle(FF), 'ob');hold on;plot(yF, angle(FT), '.r');title('FFT Phase Vs Analytic Solution');
subplot(413);plot(yF, real(FF)/scF, 'ob');hold on;plot(yF, real(FT)/anAbs, '.r');title('FFT Real Part Vs Analytic Solution');
subplot(414);plot(yF, imag(FF)/scF, 'ob');hold on; plot(yF, imag(FT)/anAbs, '.r');title('FFT Imaginary Part Vs Analytic Solution');


MSE_GSM = (sum((abs((r/max(abs(r))))-abs(an(1,:)/max(abs(an(1,:))))).^2)/sum(abs(an(1,:)/max(abs(an(1,:)))).^2))*100
MSE_AGSM = (sum((abs((s/max(abs(s))))-abs(an(2,:)/max(abs(an(2,:))))).^2)/sum(abs(an(2,:)/max(abs(an(2,:)))).^2))*100
MSE_DM = (sum((abs((dm/max(abs(dm))))-abs(an(3,:)/max(abs(an(3,:))))).^2)/sum(abs(an(3,:)/max(abs(an(3,:)))).^2))*100
MSE_FFT = (sum((abs((FF/max(abs(FF))))-abs(FT/max(abs(FT)))).^2)/sum(abs(FT/max(abs(FT))).^2))*100

% MSE_GSM = immse(r/max(abs(r)),an/max(abs(an)))
% MSE_AGSM = immse(s/max(abs(s)),an/max(abs(an)))

%% Plots for Document

% figure; 
% scs = max(abs(s));
% subplot(411);plot(y, abs(s)/scs);title('AGSM Absolute Value');
% subplot(412);plot(y, angle(s));title('AGSM Phase');
% scr = max(abs(r));
% subplot(413);plot(y, abs(r)/scr);title('GSM Absolute Value');
% subplot(414);plot(y, angle(r));title('GSM Phase');
% 
% figure;
% anAbs = max(abs(an));
% subplot(311);plot(y, abs(an)/anAbs);title('Analytic Solution Absolute Value'); 
% scr = max(abs(r));
% subplot(312);plot(y, abs(r)/scr);title('GSM Absolute Value');
% scs = max(abs(s));
% subplot(313);plot(y, abs(s)/scs);title('AGSM Absolute Value');
% 
% figure;
% subplot(311);plot(y, angle(an));title('Analytic Solution Phase'); 
% subplot(312);plot(y, angle(r));title('GSM Absolute Phase');
% subplot(313);plot(y, angle(s));title('AGSM Absolute Phase');

