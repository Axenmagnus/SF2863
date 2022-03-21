clear all; close all;clc

v = 18;
v1 = 10;
v2 = 15;
L1 = 14;
L2 = 20; mu1 = 10; mu2 = 8; h = 0.001;

PI=[];
T = 20000;
tic
for j=1:2 
PI_old=PI;
PI=[];
for i=0:3
n1=i;
n2=3-i;
T0 = 0;
T11 = 0;
T12 = 0;
T2 = 0;
t = 0;
state=0;
while t < T
    if state == 0
        t11 = exprnd(1/L1);
        t12 = exprnd(1/L2);

        if t11 < t12
            state = 11;
            t = t + t11;
            T0 = T0 + t11;
        else
            state = 12;
            t = t + t12;
            T0 = T0 + t12;
        end

    elseif state == 11;
        t0 = exprnd(1/(3 * mu1));
        t2 = exprnd(1/(L2));

        if t0 < t2
            state = 0;
            t = t + t0;
            T11 = T11 + t0;
        else
            state = 2;
            t = t + t2;
            T11 = T11 + t2;
        end

    elseif state == 12;
        t0 = exprnd(1/(3 * mu2));
        t2 = exprnd(1/L1);

        if t0 < t2;
            state = 0;
            t = t + t0;
            T12 = T12 + t0;
        else
            state = 2;
            t = t + t2;
            T12 = T12 + t2;
        end

    elseif state == 2;
        if n2==0
        t12 = exprnd(1/(3 * mu1));
        state = 12;
        t = t + t12;
        T2 = T2 + t12;
            
        elseif n1==0
        t11 = exprnd(1/(3 * mu2));
        state = 11;
        t = t + t11;
        T2 = T2 + t11;
        
        else
        t12 = exprnd(1/(n1 * mu1));
        t11 = exprnd(1/(n2 * mu2));
        
        if t11 < t12;
            state = 11;
            t = t + t11;
            T2 = T2 + t11;
        else
            state = 12;
            t = t + t12;
            T2 = T2 + t12;
        end
        end

    end
    
end
TT = T0 + T11 + T12 + T2;

pi = [T2 / TT T11 / TT T12 / TT T0 / TT];

PI= [PI; pi];
end
PI
end

Rel_error=(PI-PI_old)./PI_old
Rel_error=max(max(abs(Rel_error)))
PI=(PI+PI_old)/2

V = [0, v2, v1, v]';
Vave = PI * V
toc
