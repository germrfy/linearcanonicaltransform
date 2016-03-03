clear all;close all;
nIt = 30;
Lx = 6; N = 256; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2).*rectangularPulse(x/6);
alpha = 0.2;
A = cos(alpha*pi/2);
B = sin(alpha*pi/2);
C = -sin(alpha*pi/2);
D = cos(alpha*pi/2);
% A=1;B=1;C=2;D=3;
% A=3;B=2;C=1;D=1;
%A=0.9;B=1;C=-0.1;D=1;
%A = 0.6; B = 2; C = -0.17; D = 1.1;

yg = abs(B)/(Tx*length(g))*((0:length(g)-1) - length(g)/2);

MSE_GSM = [zeros(1,nIt)];
MSE_AGSM = [zeros(1,nIt)];
MSE_DM = [zeros(1,nIt)];
nSamp = zeros(1, nIt);
mo(nIt)= struct('cdata',[],'colormap',[]);

for i=1:nIt   
    newSize = N*(1+(0.1*(i-1)/2));
%     y = abs(B)/(Tx*newSize)*((0:newSize-1) - newSize/2);
%     f = interp1(yg,g,y, 'linear', 'extrap');
     f = padarray(g, [0, round(0.1*(i-1)*length(g))]);
     y = abs(B)/(Tx*length(f))*((0:length(f)-1) - length(f)/2);

    ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
    fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
    sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
    den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

    an = (ExAn.*(fb + sb))./(den);
    an(isnan(an)) = 0 ;

    r = gsm1d(f, Tx, A, B, C);

    s = agsm1d(f, Tx, B, C, D);
    
    d = direct_method1d(f, Tx, A, B, D);
    
%      MSE_GSM(i) = immse(r/max(abs(r)),an/max(abs(an)))*1000;
%      MSE_AGSM(i) = immse(s/max(abs(s)),an/max(abs(an)))*1000;
     MSE_GSM(i) = sum(abs((r/max(abs(r))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
     MSE_AGSM(i) = sum(abs((s/max(abs(s))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
     MSE_DM(i) = sum(abs((d/max(abs(d))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
    nSamp(i) = length(y);
    
    figure;
    sAbs = max(abs(s));
    plot(y, abs(s)/sAbs, 'b');
    hold on
    rAbs = max(abs(r));
    plot(y, abs(r)/rAbs, 'r');
    dAbs = max(abs(d));
    plot(y, abs(d)/dAbs, 'c');
    anAbs = max(abs(an));
    plot(y, abs(an)/anAbs, 'g');title('Analytic Solution Vs GSM');
    mo(i) = getframe(gcf);
    close;
end


figure;
p1 = plot([nSamp/N], MSE_GSM, '.r');
hold on
p2 = plot([nSamp/N], MSE_AGSM, 'ob');
p3 = plot([nSamp/N], MSE_DM, '-g');
legend([p1, p2, p3], 'GSM', 'AGSM', 'DM')
title('Comparison of MSE for GSM, AGSM and DM Algorithms')
xlabel('Zero Padding Factor ($${N_L \over N_0}$$)', 'Interpreter', 'Latex')
ylabel('MSE (%)')

fig = figure;
movie(fig, mo,2);


