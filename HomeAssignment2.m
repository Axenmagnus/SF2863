close all
clear
clc

lambda = [50, 40, 45, 50, 25, 48, 60, 35, 15] /1000;
T1 = 6; T2 = 8; T3 = 14; T4 = 25; T5 = 12; T6 = 18; T7 = 33; T8 = 8; T9 = 12;
T=[T1 T2 T3 T4 T5 T6 T7 T8 T9];
c1 = 12; c2 = 14; c3 = 21; c4 = 20; c5 = 11; c6 = 45; c7 = 75; c8 = 30; c9 = 22;
c=[c1 c2 c3 c4 c5 c6 c7 c8 c9];

% Initial values
%lambda = [50 40 45 51 25 48 60 35 15]/1000;
%c = [14 19 25 15 10 45 80 33 30];   
%T = [4 7 14 5 10 18 24 8 12];
C_Budget = 500;
n=length(c);

% Reasonable number of rows in R matrix
Numberofrows=8;

% Function returns matrices
[Tableau,R_matrix] = Matrices(T,lambda,n,c,Numberofrows);

% Part 1: Marginal Allocation
[Cost_vector,EBO_vector_MA,Allocation_MA,Marginal_decrease] = MarginalAllocation(T,lambda,n,c,C_Budget,Tableau);

% Part 2: Dynamic Programming
[EBO_vector_DP,Allocation_DP] = DynamicProgramming(T,lambda,n,c,C_Budget,R_matrix);

% Display Results using Marginal Allocation
disp('---------------------------------------------------------')
disp('Marginal Allocation')
disp('Budget:')
disp(C_Budget);
disp('Allocation: ')
disp(Allocation_MA(1:end,1:end)) 
disp('Expected Back-order: ')
disp(EBO_vector_MA(1,end))
disp('Total cost: ')
disp(Cost_vector(1,end))
disp('Marginal decrease: ')
disp(Marginal_decrease)

% Budgets we want to determine optimal using Dynamic Programming
C_Budgets = [0 100 150 350 500];
for i = 1:length(C_Budgets)
    [EBO,Allocation,C] = DynamicProgramming(T,lambda,n,c,C_Budgets(i),R_matrix);
    disp('---------------------------------------------------------')
    disp('Dynamic Programming')
    disp('Budget:')
    disp(C_Budgets(i));
    disp('Allocation: ')
    disp(Allocation(end,1:n)) 
    disp('Expected Back-order: ')
    disp(EBO(end,1))
    disp('Total cost: ')
    disp(C)
    
end

% Plot of Efficient curve and points
figure(1) 
plot(Cost_vector,EBO_vector_MA,'b')
hold on
plot(Cost_vector,EBO_vector_MA,'o')
legend('Efficient curve,','Efficient points') 
xlabel('Cost') 
ylabel('EBO')
xlim([0 C_Budget])
title('Marginal Allocation')

% Plot of EBO to Number of LRU purchased
figure(2)
plot(0:sum(Allocation_MA),EBO_vector_MA,'g')
xlabel('Number of LRU purchased')
ylabel('EBO')
xlim([0 sum(Allocation_MA)])
title('EBO/LRU Graph')

% Plot of marginal decrease in EBO for each new unit of cost spent
figure(3)
plot(Cost_vector(1,1:end-1), -Marginal_decrease,'m')
xlabel('Cost')
ylabel('Marginal decrease in EBO')
xlim([0 Cost_vector(1,end-1)])
title('Marginal decrease in EBO for each new unit of cost spent')

% Plot of dynamic programming solutions
figure(4) 
plot(0:C_Budget,EBO_vector_DP,'r')
xlabel('Cost') 
ylabel('EBO')
xlim([0 C_Budget])
legend('Dynamic Progamming solutions')
title('Dynamic Progamming')

% Plot of MA and DP solutions for comparison
figure(5)
plot(Cost_vector,EBO_vector_MA,'b') 
hold on
plot(0:C_Budget,EBO_vector_DP,'r')
hold on 
plot(Cost_vector,EBO_vector_MA,'o')
xlabel('Cost')
ylabel('EBO')
xlim([0 C_Budget])
legend('Marginal Allocation','Dynamic Programming')
title('Comparison')


% Computing R matrix and Tableau with R(i,j)/c(j)
function [Tableau,R_Matrix] = Matrices(T,lambda,n,c,Numberofrows)
    R = ones(1,n);
    for i=1:Numberofrows
        Rx = R(end,:) - (lambda.*T).^(i-1)./factorial(i-1).*exp(-lambda.*T);  %/factorial(s-1)*exp(-landa*T)
        R = [R;Rx];
    end
    
    R_Matrix = R(2:end,:);
    Tableau = R_Matrix./c;
end    

% Marginal Allocation
function [Cost_vector,EBO_vector,s,Marginal_EBO] = MarginalAllocation(T,lambda,n,c,C_Budget,Tab)
    s = zeros(1,n); % Initial allocation
    EBO = lambda.*T; % Initial EBO
    C = 0; % Initial cost
    
    % Initiating vectors
    EBO_vector = sum(EBO);        
    Cost_vector = C;
    
    % The marginal decrease in EBO for each new unit of cost
    Marginal_EBO=[];
    

    % Infinity loop with break if we can't afford anything
    while true
        % Check the cost effectiveness for
        % every LRU for next purchase
        Cost_effectiveness = [];
        for j = 1:n
            index = s(j)+1;
            Cost_effectiveness = [Cost_effectiveness, Tab(index,j)]; 
        end
        
        % Affordable check
        for i=1:n
            if c(i)>C_Budget-C
                Cost_effectiveness(1,i)=0;
            end 
        end
        
        % Break loop if we can't afford anything
        if sum(Cost_effectiveness,'all') == 0
            break
        end
        

        if C_Budget<C
            break
        end
        % Searching index for most cost-effective purchase
        m = find(Cost_effectiveness==max(Cost_effectiveness));
        
        % If we find two or more, we choose one
        if length(m)> 1
            m=m(1);   
        end
        
        % Updating LRU allocation, EBO, marginal decrease and cost
        s(m) = s(m)+1;
        EBO(m)=EBO(m) - max(Cost_effectiveness)*c(m);
        EBO_vector = [EBO_vector,sum(EBO)];
        Marginal_EBO = [Marginal_EBO,Cost_effectiveness(m)];
        C = s*c';
        Cost_vector = [Cost_vector, C];
    end
    
end


% Dynamic Programming
function [EBO_vector,s,optimal_cost] = DynamicProgramming(T,lambda,n,c,C_Budget,R)
    s = zeros(1,n); % Initial allocation
    s_update = zeros(1,n); % Allocation update vector
    cost_stage = 1; % Initial cost stage
    EBO_vector = sum(lambda.*T);
    
    % Infinity loop with break if cost budget is reached
    while true
        if C_Budget == 0
            optimal_cost=0;
            break
        end
        
        % n+1 decisions, purchasing a LRU or not purchasing 
        delta_EBO = inf(1,n+1);
        
        % Not purchasing => EBO remains
        delta_EBO(1) = EBO_vector(cost_stage);
        
        % Checking if we can afford LRU(i), then computing
        % the optimal allocation if we purchase LRU(i)
        for i = 1:n
            if cost_stage-c(i) >= 0 
                s1 = s(cost_stage-c(i)+1,:);
                if c*s1'<= cost_stage - c(i)       
                    delta_EBO(i+1) =EBO_vector(cost_stage-c(i)+1)-R(s1(i)+1,i);
                end
            end
        end
        
        % Finding the minimal EBO allocation of the optimal allocations
        minimal_EBO = min(delta_EBO);
        
        x = find(minimal_EBO == delta_EBO,1);
        
        % If we purchase, we update the allocation
        if x ~= 1
            s_update = s(cost_stage-c(x-1)+1,:);
            s_update(x-1) = s_update(x-1)+1; 
        end
        
        % EBO vector and allocation update
        EBO_vector = [EBO_vector;minimal_EBO];
        s = [s;s_update];
        
        cost_stage = cost_stage+1;
        
        if cost_stage == C_Budget+1
            optimal_cost=s(end,1:n)*c';
            break  
        end
    end
end

