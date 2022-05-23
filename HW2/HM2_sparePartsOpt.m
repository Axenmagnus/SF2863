clc, close all, clear all;

lambda = [50, 40, 45, 50, 25, 48, 60, 35, 15] /1000;
T1 = 6; T2 = 8; T3 = 14; T4 = 25; T5 = 12; T6 = 18; T7 = 33; T8 = 8; T9 = 12;
T=[T1 T2 T3 T4 T5 T6 T7 T8 T9];
c1 = 12; c2 = 14; c3 = 21; c4 = 20; c5 = 11; c6 = 45; c7 = 75; c8 = 30; c9 = 22;
c=[c1 c2 c3 c4 c5 c6 c7 c8 c9];
Cmax=500;

MA=[];

for i=1:9
   MA=[MA R(lambda,T,0,i)/c(i)];
end

ind=zeros(1,9);

C=0;
EBO=dot(lambda,T);
C_vek=[];
EBO_vek=[];
dEBO=[];

while C<=Cmax
    C_vek=[C_vek C];
    EBO_vek= [EBO_vek EBO];
    [M l]=max(MA);
    if M==0
        break
    end
    
    C=C+c(l);
    while C>=Cmax %Cant add the requested spare part
        MA(l)=0;
        C=C-c(l);
        [M l]=max(MA);
        C=C+c(l);
        if M==0
            break
        end
    
    end
        EBO=EBO-R(lambda,T,ind(l),l);

        ind(l)=ind(l)+1;
        MA(l)= R(lambda,T,ind(l),l)/c(l);
        
    
end


figure(1)
plot(C_vek,EBO_vek);
hold on
plot(C_vek,EBO_vek,'*');
xlabel('Cost')
ylabel('EBO')
legend('Efficient Curve','Efficient Points' )
grid on
hold on
saveas(figure(1), [pwd append('/curve.png')])


dEBO
max(ind);

MA=[];
for j=0:max(ind)
   MA_row=[]; 
for i=1:9
  MA_row=[MA_row R(lambda,T,j,i)/c(i)];
end
MA=[MA; MA_row];
end
MA



%% Dynamic programming
n=length(lambda);
Cmax=500;
optimalsolutions=zeros(Cmax,n);
EboMatrix=10*ones(Cmax,1);
EboMatrix(1)=sum(lambda.*T);


for cStage=1:Cmax
    if cStage<min(c)    
        EboMatrix(cStage)=sum(lambda.*T);
    else
        TempEbo=zeros(n,1);
        for i=1:n
            tempcost=cStage-c(i);
            if tempcost>=0
                PreviousSolution=optimalsolutions(tempcost+1,:);
                Totalcost=sum(c.*PreviousSolution);
                %Totalcost
                if Totalcost<=tempcost;
                    TempEbo(i)=EboMatrix(tempcost+1,:)-R(lambda,T,optimalsolutions(tempcost+1,i),i);
                end

            end
        end

        for i=1:length(TempEbo)
            if TempEbo(i)==0
                TempEbo(i)=inf; 
            end
        end 
        
        [MinEbo,l]=min(TempEbo);
        if MinEbo<min(EboMatrix)
       
        AddedVal=zeros(n,1)';
        AddedVal(l)=1;
        if cStage-c(l)>=0
            optimalsolutions(cStage,:)=optimalsolutions(cStage-c(l)+1,:)+AddedVal;
            EboMatrix(cStage)=TempEbo(l);

        end
        else
            EboMatrix(cStage)=EboMatrix(cStage-1);
            optimalsolutions(cStage,:)=optimalsolutions(cStage-1,:);
        end
    end
    

end
EboMatrix=[sum(lambda.*T); EboMatrix];
sum(optimalsolutions(500,:).*c)
hold on
plot(0:Cmax, EboMatrix)
legend('Efficient Curve','Efficient Points','Dynamic Programming')
saveas(figure(1), [pwd append('/Dcurve.png')])


%%
function R=R(lambda,T,i,j)
sum=0;
for k=0:i  
    sum=sum+((lambda(j)*T(j))^k * exp(-lambda(j)*T(j)))/factorial(k);
end
R=1-sum;
end
