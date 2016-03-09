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
MSE_DM = [zeros(1,nIt)];
nSamp = zeros(1, nIt);
mo(nIt)= struct('cdata',[],'colormap',[]);

for i=1:nIt     
    f = padarray(g, [0, round(0.1*(i-1)*length(g))]);
    yGSM = abs(A)*Tx*((0:length(f)-1) - length(f)/2);       %output vector GSM
    yAGSM = (Tx/abs(D))*((0:length(f)-1) - length(f)/2);    %output vector AGSM
    y = [yGSM; yAGSM];

    ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
    fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
    sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
    den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

    an = (ExAn.*(fb + sb))./(den);
    an(isnan(an)) = 0 ;

    r = gsm1d(f, Tx, A, B, C);

    s = agsm1d(f, Tx, B, C, D);
        
    MSE_GSM(i) = sum(abs((r/max(abs(r))-(an(1,:)/max(abs(an(1,:)))))).^2)/sum(abs(an(1,:)/max(abs(an(1,:)))).^2)*100;
    MSE_AGSM(i) = sum(abs((s/max(abs(s))-(an(2,:)/max(abs(an(2,:)))))).^2)/sum(abs(an(2,:)/max(abs(an(2,:)))).^2)*100;
%      MSE_GSM(i) = sum(abs((r/max(abs(r))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
%      MSE_AGSM(i) = sum(abs((s/max(abs(s))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
    nSamp(i) = length(y);
    
%     figure;
%     rAbs = max(abs(r));
%     plot(y(1,:), abs(r)/rAbs, 'b');
%     hold on
% %     rAbs = max(abs(r));
% %     plot(y(2,:), abs(r)/rAbs, 'r');
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
title('Comparison of MSE for GSM, AGSM Algorithms')
xlabel('Zero Padding Factor ($${N_L \over N_0}$$)', 'Interpreter', 'Latex')
ylabel('MSE (%)')
% 
% figure;
% axis off;
% movie(mo, 12);

