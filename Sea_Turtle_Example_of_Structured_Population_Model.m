%% Sea_Turtle_Example_of_Structured_Population_Model

% We should first construct the projection matrix A
% Before that, we have make two assumptions: 
% 1st: All the parameters in A do not change over time, at least during the
% observation
% 2st: All parameters is uniform within each stage among all ages
% A(i,i) = 1-1/(1+pi+...pi^(di-1))
% A(i+1,i) = pi^(di)/(1+pi+...pi^(di-1))

% Define projection matrix A
A=[ 0     0     0     0     127      4      80  ; 
   .6747  .737  0     0       0      0       0  ; 
    0  .0486 .6610   0       0      0       0  ; 
    0     0  .0147  .6907    0      0       0  ;
    0     0    0    .0518    0      0       0  ;
    0     0    0     0      .8091   0       0  ;
    0     0    0     0       0     .8091   .8089];

%Initial population size in each stage class
n0 = [900;
      900;
      160;
       10;
        5;
        4;
       20];    
                
tmax    = 100;
nt      = zeros(7,tmax);
nt(:,1) = n0 ;

for t=2:tmax
   nt(:,t)=A*nt(:,t-1) ;    
end

figure
set(gca,'FontSize',20)
plot(1:tmax,nt','.-','MarkerSize',14)
xlabel('t','FontSize',20)
ylabel('n','FontSize',20)
legend('stage 1','stage 2','stage 3','stage 4','stage 5','stage 6','stage 7')

ppp = polyfit(25:tmax,log(nt(1,25:tmax)),1);
lambda_estimate=exp(ppp(1));

% sensitivity analysis
% Sij = partial lambda/partial aij = vi*wj/sum(vk.*wk)
% where v is the eigenvector of A' associated with lambda and w is the
% eigenvector of A associated with lambada
% In short S = v*w'/sum(v.*w)
% While sensitivity reflects the absolute change, % elasticity reflects 
% the percentage change
% Eij = (partial lambda/lambda) / (partial aij/aij) 
% = partial log(lambda)/ partial log(aij)
% = aij/lambda*Sij
% In short E = A.*S/lambda

% Compute right eigenvalues and right eigenvectors(the default ones) of A
[W,D_r]  = eig(A);
% Compute left eigenvalues and left eigenvectors of A
[V,D_l] = eig(A');
 
dr = diag(D_r);
[lambda,indr] = max(dr); % find dominant eigenvalue
w  = W(:,indr); % corresponding (dominant) eigenvector

dl = diag(D_l);
[lambda,indl] = max(dl); % find dominant eigenvalue
v  = V(:,indl); % corresponding eigenvector

S = v*w'/sum(v.*w); % snesitivity
E = A.*S/lambda; % elasticity

% elasticities for P, F and G
for i=1:6
    elastg(i)=E(i+1,i);
end
elastp=diag(E);
elastf=E(1,:);

figure;
plot(elastp,'^-','MarkerSize',8,'LineWidth',1);
hold on;
plot(elastg,'*-','MarkerSize',8,'LineWidth',1);
hold on;
plot(elastf,'o-','MarkerSize',8,'LineWidth',1)
legend('P_i','G_i','F_i');
xlabel('Stage');
ylabel('Elasticity');