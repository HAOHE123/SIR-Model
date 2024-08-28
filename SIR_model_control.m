clearvars;
close all

N=8; %how many nodes in the network
K=3; %how many nodes you want to control %aftermytest,this should be NbyN
total=10;%how many runs

%Error using dare (line 101)
%Unable to solve the specified Riccati equation because the Symplectic spectrum is too close
%to the unit circle.
%we use try catch to skip it
failed_run=0;


time=20;
L_total=zeros(N,total);
E=zeros(total,time);


x = zeros(N,time);
x(:,1)= ones(N,1);
x_desired=zeros(N,1);


%If we want to make sure how many nodes need to be controlled, we need
%to make that row become all 0
%m = randi(1,n)
%v = randi([m,1],1,n)
%mv cannot be used


%m should be randi(N,N)
%v should be randi([1,N],1,N) still not work since it must be larger than 1

for total_run=1:total
e=zeros(1,time);
    %% Random A and B
A=eye(N);
for i=1:N
    for j=1:N
        k=rand;
        if k>0.5
            k=1;
        else
            k=0;
        end
        A(i,j)=double(k);
    end
end

A=A-diag(diag(A));
A=triu(A)+(triu(A))';
% A=A+eye(N); % THIS WAS WRONG GIVEN THE DEFINITION OF THE TRANS_FUNCTION
B=zeros(N,N);
node_control_count=zeros(N,1);

for i=1:K
    
    position=round(rand*N);  %random position and check if contains same node
    for l=1:1000
        if any(node_control_count==position) 
            position=round(rand*N);
        else
            node_control_count(i)=position;
            break;
        end
    end
     k=rand;
     if k>0.5 % THIS PART IS A BIT POINTLESS GIVEN THAT THE GAIN WILL JUST PUT - or + when necessary
        k=1;
     else
        k=-1;
     end
        B(i,position)=double(k);
end


%% dare function with Q and R

R=eye(N);
Q=10*R;
try
[X,L,G]=dare(A,B,Q,R);
catch 
    failed_run=failed_run+1;
    continue;
end


%% Error Estimation


for t=1:time-1
    trans_function{t} = double(logical(A^t+A^(t-1)));
     
    x(:,t+1) = A*x(:,t)-B*G*x(:,t);
    
    e(t)=norm(x(:,t)-x_desired); % IT IS BETTER TO USE A NORM THAN A SUBTRACTION TO AVOID LATER CANCELATION
    
end

e(time)=norm(x(:,time)-x_desired); % WE ALSO WANT TO COMPUTE THE ERROR OF THE FINAL STEP

L_total(:,(total_run-failed_run))=L;
E(total_run,:)=e;
end

% NO NEED TO DO THIS PART OF THE CODE
% %% L average & e average
% %L average
% L_ave=zeros(N,1);
% for i=1:(total_run-failed_run)
%     L_ave=L_ave+L_total((total_run-failed_run));
% end
% L_ave=L_ave/(total_run-failed_run);
% 
% %e average
% e_ave=zeros(N,time);
% for i=1:(total_run-failed_run)
%     for j=1:time
%         e_ave(:,j)=E(:,j+time*((total_run-failed_run)-1));
%     end
% end
% e_ave=e_ave/(total_run-failed_run);


%% Draw Figures

figure(1);
boxplot(sort(L_total)')
title('Statistics of the eigenvalues of the controlled system');
xlabel('$i$','Interpreter','latex');
ylabel('$\lambda_i$','Interpreter','latex');
saveas(1,['eigenvalues-statistics-N-' num2str(N) '-C-' num2str(K)],'pdf');    

figure(2);
boxplot(E);
title('Statistics of the error over time');
xlabel('$k$','Interpreter','latex');
ylabel('$\|x(k)\|$','Interpreter','latex');
saveas(2,['error-statistics-N-' num2str(N) '-C-' num2str(K)],'pdf');    
