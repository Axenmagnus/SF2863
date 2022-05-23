clear all;close all; clc
Solutions = [];

v = 18;
v1 = 10;
v2 = 15;
L1 = 14;
L2 = 20;
mu1 = 10; mu2 = 8;
h = 0.001;
Q = @(n1, n2) [- (n1 * mu1 + n2 * mu2), n2 * mu2, n1 * mu1, 0;
        L2, - (L2 + 3 * mu1), 0, 3 * mu1;
        L1, 0, - (L1 + 3 * mu2), 3 * mu2;
        0, L1, L2, - (L1 + L2)];

U = zeros(4, 1);
U(end) = 1;

for i = 0:3
    n1 = i;
    n2 = 3 - i;
    c = Q(n1, n2);
    c(:, end) = 1;
    Sol = c' \ U;
    Solutions = [Solutions; Sol'];
        
end

A=Q(n1, n2); % doesnt matter which strategy
A=A(2:end,2:end);
b=-[1/(-A(1,1)); 1/(-A(2,2)); 1/(-A(3,3))];
A=[A(1,:)/(-A(1,1)); A(2,:)/(-A(2,2)); A(3,:)/(-A(3,3))];
t=A\b 
Solutions

SUM = [];

for j = 1:length(Solutions')
    sum = 0;

    for i = 1:length(Solutions')
        sum = sum + Solutions(j, i);
    end

    SUM = [SUM sum];
end

SUM
L = Q(3, 0);
Solutionss = L \ U;

V = [0, v2, v1, v]';

Vave = Solutions * V;
