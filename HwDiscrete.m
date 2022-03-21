clear all; close all;clc

v = 18;
v1 = 10;
v2 = 15;
L1 = 14;
L2 = 20; mu1 = 10; mu2 = 8; h = 0.001;
Q = @(n1, n2) [- (n1 * mu1 + n2 * mu2), n2 * mu2, n1 * mu1, 0;
        L2, - (L2 + 3 * mu1), 0, 3 * mu1;
        L1, 0, - (L1 + 3 * mu2), 3 * mu2;
        0, L1, L2, - (L1 + L2)];


N=10*10^6;
PI=[];
cond=1;
format short
% To gain averages
while cond==1
    tic
for j=1:2
    
    PI_old=PI;
    PI=[];
for i=0:3
n1=i;
n2=3-i;

P = Q(n1, n2) * h + eye(length(Q(n1, n2)));

t=0;
sum0=0;
sum1=0;
sum2=0;
sum3=0;
state=0;
States=[];
for k=0:N
    val=rand;
    %States=[States,state];
    
    if state==0
        %Consider First row in pmatrix
        
        if val<P(1,1)
            state=0;
            sum0=sum0+1;
            t=k*h;
        elseif P(1,1)<val && val<(P(1,1)+P(1,2))
            state=1;
            sum0=sum0+1;
            t=k*h;
        elseif (P(1,1)+P(1,2))<val && val<(P(1,1)+P(1,2)+P(1,3))
            state=2;
            sum0=sum0+1;
            t=k*h;
        end
    elseif state==1
        if val<P(2,1)
            state=0;
            sum1=sum1+1;
            t=k*h;
        elseif P(2,1)<val && val<P(2,1)+P(2,2)
            state=1;
            sum1=sum1+1;
            t=k*h;
        elseif P(2,1)+P(2,2)<val && val<P(2,1)+P(2,2)+P(2,4)
            state=3;
            sum1=sum1+1;
            t=k*h;
        end
    elseif state==2
        if val<P(3,1)
            state=0;
            sum2=sum2+1;
            t=k*h;
        elseif P(3,1)<val && val<P(3,1)+P(3,3)
            state=2;
            sum2=sum2+1;
            t=k*h;
        elseif P(3,1)+P(3,3)<val && val<P(3,1)+P(3,3)+P(3,4)
            state=3;
            sum2=sum2+1;
            t=k*h;
        end
    elseif state==3
        if val<P(4,2)
            sum3=sum3+1;
            state=1;
            t=k*h;
        elseif P(4,2)<val && val<P(4,2)+P(4,3)
            sum3=sum3+1;
            state=2;
            t=k*h;
        elseif P(4,2)+P(4,3)<val && val<P(4,2)+P(4,3)+P(4,4)
            sum3=sum3+1;
            state=3;
            t=k*h;
        end
    end
    
            
    
    
    
end

pi0=sum0/N;
pi1=sum1/N;
pi2=sum2/N;
pi3=sum3/N;

pi=[pi0 pi1 pi2 pi3];
PI=[PI; pi];
end


end

Rel_error=(PI-PI_old)./PI_old;
Rel_error=max(max(abs(Rel_error)))
PI=(PI+PI_old)/2;

V = [0, v2, v1, v]';
Vave = PI * V
if Rel_error<0.01
    cond=2;
else
    cond=1;
end
N=N+1000000
toc
end

