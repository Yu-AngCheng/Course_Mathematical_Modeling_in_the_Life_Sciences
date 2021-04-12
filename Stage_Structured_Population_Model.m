%% Stage-Structured Population Model
clear all
close all

% n(a,t) number of females at age a & time t
% Amax maximum age
% tmax maximum time
% p(a) probability of age a surviving to age a+1
% f(a) number of new-borns per year per age-a female

% A Leslie maxtrix(projection matrix),where aij means the stage j at time t
% transfers to stage i at time t+1

% dynamic equation
% n(:,t+1) = A*n(:,t)
% if A has single real positive dominant eigenvalue, then n(:,t) =
% c*lambda^(t)*w where w is the corresponding eigenvector

% if A is k*k, then A is power positive if and only if A^(k^2-2k+2) is
% positive

% Perron-Frobenius Theorem:
% if A is non-negative, square, power positive, then A has single real
% positive dominant eigenvalue, and the corresponding eigenvector is
% also real and positive

% data here
Amax = 2;
tmax = 100;
p = [0.5,0.25,0];
f = [0,1,5];
n = NaN(Amax,tmax);
n(1,1) = 100; n(2,1) = 50; n(3,1) = 13;
A = [0,  1,  5;
     .5,  0,  0;
     0, .25, 0];
 
 for t = 2:tmax
     n(:,t) = A*n(:,t-1);
 end
 
 % plot to show
figure;
for a = 1:Amax+1
    plot(1:tmax,n(a,:),'LineWidth',1);
    hold on;
end
legend({'Age 0','Age 1','Age 2'});

figure;
for a = 1:Amax+1
    semilogy(1:tmax,n(a,:),'LineWidth',1);
    hold on;
end
legend({'Age 0','Age 1','Age 2'});

