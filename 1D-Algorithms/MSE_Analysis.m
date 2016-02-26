clear all;
nIt = 20;
Lx = 6; N = 256; Tx = Lx/N; fs = 1/Tx;
x = [-3:Tx:3 - Tx];
g = exp(-x.^2).*rectangularPulse(x/6);
alpha = 0.3;
A = cos(alpha*pi/2);
B = sin(alpha*pi/2);
C = -sin(alpha*pi/2);
D = cos(alpha*pi/2);
% A=0.9;B=-1;C=0.1;D=1;
% A=3;B=2;C=1;D=1;

MSE_GSM = [zeros(1,nIt)];
MSE_AGSM = [zeros(1,nIt)];
nSamp = zeros(1, nIt);

for i=1:nIt
    f = padarray(g, [0, round(0.1*(i-1)*length(g))]);
    y = abs(B)/(Tx*length(f))*((0:length(f)-1) - length(f)/2);

    ExAn = -1i*sqrt(1i/B)*exp((pi*(1i*B*D+(-1+A*D)*pi)*y.^2)/(B*(B-1i*A*pi)))*sqrt(pi);
    fb = (3*B-3*1i*A*pi+1i*pi*y).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))).*erfz(sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))));
    sb = sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*(3*B-1i*pi*(3*A+y)).*erfz(sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi))));
    den = 2*(B-1i*A*pi)*sqrt(((3*B-3*1i*A*pi+1i*pi*y).^2)/(B*(B-1i*A*pi))).*sqrt(((3*B-1i*pi*(3*A+y)).^2)/(B*(B-1i*A*pi)));

    an = (ExAn.*(fb + sb))./(den);

    r = gsm1d(f, Tx, A, B, C);

    s = agsm1d(f, Tx, B, C, D);
    
%     MSE_GSM(i) = immse(r/max(abs(r)),an/max(abs(an)));
%     MSE_AGSM(i) = immse(s/max(abs(s)),an/max(abs(an)));
    MSE_GSM(i) = sum(abs((r/max(abs(r))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
    MSE_AGSM(i) = sum(abs((s/max(abs(s))-(an/max(abs(an))))).^2)/sum(abs(an/max(abs(an))).^2)*100;
    nSamp(i) = length(y);
end


figure;
p1 = plot([nSamp/N], MSE_GSM, 'or');
hold on
p2 = plot([nSamp/N], MSE_AGSM, '.b');
legend([p1, p2], 'SM', 'ASM')
title('Comparison of MSE for GSM and AGSM')
xlabel('$${N_L \over N_0}$$', 'Interpreter', 'Latex')
ylabel('MSE')
