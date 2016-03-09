function SBRAnal(A, B, C, D, R)

i = size(R);
if i(2) > 1 
    R = R';
end

tol = eps(4);
if abs(det([A B; C D])-1) > tol
    error('LCT Error - Determinant constraint not met')
end

X = [abs((-A*sqrt(R)+B*sqrt(1./R))).*abs(-C*sqrt(R)+D*sqrt(1./R)) abs((A*sqrt(R)+B*sqrt(1./R))).*abs(C*sqrt(R)+D*sqrt(1./R)) abs((A*sqrt(R)+B*sqrt(1./R))).*abs(-C*sqrt(R)+D*sqrt(1./R)) abs((-A*sqrt(R)+B*sqrt(1./R))).*abs(C*sqrt(R)+D*sqrt(1./R))];
d = [abs(1+(A/B).*R) abs(1-(A/B).*R) X]
s = [abs(1+(B/A)*(1./R)) abs(1-(B/A)*(1./R)) X]
as = [abs(1+(C/D)*R) abs(1-(C/D)*R) X]
DM = max([abs(1+(A/B).*R) abs(1-(A/B).*R) X], [], 2)
SM = max([abs(1+(B/A)*(1./R)) abs(1-(B/A)*(1./R)) X], [], 2)
ASM = max([abs(1+(C/D)*R) abs(1-(C/D)*R) X], [], 2)

 p1 = plot(R, DM, '.r');
 
 hold on
 
 p2 = plot(R, SM, '*b');
 
 p3 = plot(R, ASM, 'og');
 
 legend([p1, p2, p3], 'DM', 'GSM', 'AGSM');
 xlabel('Space Bandwidth Ratio')
 ylabel('$${N_L \over N_0}$$', 'Interpreter', 'Latex')
 