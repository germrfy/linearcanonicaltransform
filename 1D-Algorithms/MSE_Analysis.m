clear all;close all;
nIt = 15;
Lx = 6; Ly = 6; N = 32; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2).*rectangularPulse(x/6);
% alpha = 0.2;
% A = cos(alpha*pi/2);
% B = sin(alpha*pi/2);
% C = -sin(alpha*pi/2);
% D = cos(alpha*pi/2);
% A=1;B=1;C=2;D=3;
% A=3;B=2;C=1;D=1;
% A=0.9;B=1;C=-0.1;D=1;     %%PAPER
%A = 0.6; B = 2; C = -0.17; D = 1.1;
A=0.5;B=4;C=-0.125;D=1;     %%PAPER

yg = abs(B)/(Tx*length(g))*((0:length(g)-1) - length(g)/2);

MSE_GSM = [zeros(1,nIt)];
MSE_AGSM = [zeros(1,nIt)];
MSE_FFT = [zeros(1,nIt)];
nSamp = zeros(1, nIt);
mo(nIt)= struct('cdata',[],'colormap',[]);

for i=1:nIt     
    f = padarray(g, [0, round(0.1*(i-1)*length(g))]);
    yGSM = abs(A)*Tx*((0:length(f)-1) - length(f)/2);       %output vector GSM
    yAGSM = (Tx/abs(D))*((0:length(f)-1) - length(f)/2);    %output vector AGSM    
    yF= ((2*pi*fs)/(length(f)))*((0:length(f)-1) - length(f)/2); %output for FFT
    y = [yGSM; yAGSM];                                      % AGSM and GSM vectors
                            %LCT analytic Soln
    ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
    fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
    sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
    den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

    an = (ExAn.*(fb + sb))./(den);
    an(isnan(an)) = 0 ;
                            % FT analytic Soln
    FT = (exp(-(yF.^2/4)).*(erfz(3-(1i*yF/2))+erfz(3+(1i*yF/2))))/(2*sqrt(2));

    r = gsm1d(f, Tx, A, B, C);

    s = agsm1d(f, Tx, B, C, D);
    
    FFT = fftshift(fft(fftshift(f)));
        
    MSE_GSM(i) = sum(abs((r/max(abs(r))-(an(1,:)/max(abs(an(1,:)))))).^2)/sum(abs(an(1,:)/max(abs(an(1,:)))).^2)*100;
    MSE_AGSM(i) = sum(abs((s/max(abs(s))-(an(2,:)/max(abs(an(2,:)))))).^2)/sum(abs(an(2,:)/max(abs(an(2,:)))).^2)*100;
    MSE_FFT(i) = sum(abs((FFT/max(abs(FFT))-(FT/max(abs(FT))))).^2)/sum(abs(FT/max(abs(FT))).^2)*100;
    nSamp(i) = length(y);
    
%     figure;
%     rAbs = max(abs(r));
%     plot(y(1,:), abs(r)/rAbs, 'b');
%     hold on
%     rAbs = max(abs(r));
%     plot(y(2,:), abs(r)/rAbs, 'r');
%     anAbs = max(abs(an(1,:)));
%     plot(y(1,:), abs(an(1,:))/anAbs, 'g');title('Analytic Solution Vs GSM');
%     mo(i) = getframe(gcf);
%     close;
end


figure;
p1 = plot([nSamp/N], MSE_GSM, 'r*');
hold on
p2 = plot([nSamp/N], MSE_AGSM, 'bo');
legend([p1, p2], 'GSM', 'AGSM')
title('MSE Analysis for GSM and AGSM Algorithms')
xlabel('Zero Padding Factor', 'Interpreter', 'Latex')
ylabel('MSE (%)')
% 
figure;
plot([nSamp/N], MSE_FFT, 'r*')
title('MSE Analysis for FFT Function in MATLAB')
xlabel('Zero Padding Factor', 'Interpreter', 'Latex')
ylabel('MSE (%)')
