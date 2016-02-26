clear all;
close all;
%A = 0.9; B = -1; C = 0.1; D = 1;
A = 0.6; B = 2; C = -0.17; D = 1.1;
% t = 2*pi/6;
% A = cos(t); D = A;
% B = sin(t); C = -B;
% alpha = 0.5;
% A = cos(alpha*pi/2);
% B = sin(alpha*pi/2);
% C = -sin(alpha*pi/2);
% D = cos(alpha*pi/2);

R = [0:.1:10]';

%tol = eps(1);
%abs(det([A B; C D])-1) > tol
SBRAnal(A, B, C, D, R);
%x = det([A B; C D])

%x - 1




DM = max([abs(1+(A/B).*R) abs(1-(A/B).*R)], [], 2);
SM = max([abs(1+(B/A)*(1./R)) abs(1-(B/A)*R)], [], 2);
ASM = max([abs(1+(C/D)*R) abs(1-(C/D)*R)], [], 2);