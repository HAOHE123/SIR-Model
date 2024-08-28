clear all
close all

nb_experiments = 1000;

success_runs = 0;
failed_runs = 0;

success_runs2 = 0;
failed_runs2 = 0;

nb_initial_nodes = zeros(1,nb_experiments);

for tests=1:nb_experiments
N=6;
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
% A=[0 1 1 0 0 0;
%    1 0 0 1 1 0;
%    1 0 0 0 0 1;
%    0 1 0 0 0 0;
%    0 1 0 0 0 0;
%    0 0 1 0 0 0];
%suppose the observor nodes are node 1, 2 and 3, 4
C=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0];
C=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 1 0 0];
N=size(A,1);
K=size(C,1);
time = 20;

%observors
x = zeros(N,time);
y = zeros(size(C,1),time);

%suppose the active node is 6
x(:,1)=[0;0;0;0;0;1];
O=zeros(size(C,1)*N,N);
O(1:size(C,1),:)=C;
trans_function=cell(1,time);


for t=1:time
    trans_function{t} = double(logical(A^t+A^(t-1)));
    if t <= N-1
        O(1+K*t:K*(t+1),:)=C*trans_function{t};
    end
    x(:,t+1) = trans_function{t}*x(:,1);
    y(:,t) = C*x(:,t);
end

node=zeros(N,time); 

% X0 reconstructed
Y = reshape(y(:,1:N),[size(C,1)*N,1]);
X0 = O'*O\O'*Y;

X = sdpvar(N,1);
F = [O*X == Y];
optimize(F,norm(X,1));
X=double(X)


if sum(sum(abs(x(:,1)-double(X)))) < 1e-3
    success_runs = success_runs +1;
else
    failed_runs = failed_runs +1;
end







x2 = zeros(N,time);
y2 = zeros(size(C,1),time);
x2(:,1)=[0;0;0;0;0;0];

%impulse node 6 at time 5
ut=zeros(N,time);
impulse=[0;0;0;0;0;1];
ut(:,5)=impulse;

%observers
O=zeros(size(C,1)*time,N);
O(1:size(C,1),:)=C;
trans_function=cell(1,time);
p=zeros(N,time);
phase=zeros(1,6);

for t=1:time-1
    trans_function{t} = double(logical(A^t+A^(t-1)));
    if t <= time-1
        O(1+K*t:K*(t+1),:)=C*trans_function{t};
    end
    
    for i=1:t
        %x2 is z
        phase=trans_function{t-i+1}*ut(:,i);
        p(:,t)=phase+p(:,t);
        p(find(p>0))=1;
    end 
    
 
    x2(:,t+1) = trans_function{t}*x2(:,1)+p(:,t)+ut(:,t+1);
    y2(:,t) = C*x2(:,t);
end

y2(:,time)=C*x2(:,time);
X2=zeros(N,1);

for i=1:time-1
% X2 reconstructed

if length(intersect(find(X2<1e-5),find(X2>-1e-5)))<5-1e-5 || length(find(X2>1-1e-5))<1-1e-5 || sum(double(X2))>1+1e-5
%if length(intersect(find(X2<1e-3),find(X2>-1e-3)))<5-1e-3 || length(find(X2>1-1e-3))<1-1e-3

%more if constraints, the process will be much more accurate but also much
%more slower

Y2 = reshape(y2(:,i:time),[size(C,1)*(time-i+1),1]);
O2=O(1:time*K-K*(i-1),:);
X2 = sdpvar(N,1);
F = [O2*X2 == Y2];
optimize(F,norm(X2,1));
X2=double(X2);
end
if sum(sum(abs(impulse-double(X2)))) < 1e-3
    break;
end
end

%X2 is in fact the impluse, which is the source but not u(0)
%in the new model, in fact we surly don`t need O since there is no X and we
%only have a new state z=x+.$$$.u
if sum(sum(abs(impulse-double(X2)))) < 1e-3
    success_runs2 = success_runs2 +1;
else
    failed_runs2 = failed_runs2 +1;
end

end
