clear all;
close all;
%A = 0.9; B = -1; C = 0.1; D = 1;
%A = 0.6; B = 2; C = -0.17; D = 1.1;
% t = 2*pi/6;
% A = cos(t); D = A;
% B = sin(t); C = -B;
% alpha = 0.9;
% A = cos(alpha*pi/2);
% B = sin(alpha*pi/2);
% C = -sin(alpha*pi/2);
% D = cos(alpha*pi/2);
%A=0.01;B=1;C=-0.99;D=1;
%A=0.5;B=4;C=-0.125;D=1;
%A=1;B=4;C=-0.125;D=7;

A=1;B=2;C=2.5;D=6;
A=0.1;B=0.9;C=-1;D=1;
R = [[0:0.1:1.5] [1.5:0.25:6]]';

%tol = eps(1);
%abs(det([A B; C D])-1) > tol
SBRAnal(A, B, C, D, R);
%x = det([A B; C D])

%x - 1
