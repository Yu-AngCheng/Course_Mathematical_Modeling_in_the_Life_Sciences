%% Age-Structured Population Model
clear all
close all

% n(a,t) number of females at age a & time t
% Amax maximum age
% tmax maximum time
% p(a) probability of age a surviving to age a+1
% f(a) number of new-borns per year per age-a female

% original equation
% n(a,t+1) = p(a-1)*n(a-1,t)
% n(0,t+1) = sum a (f(a).*n(a,t),1)
% simplified equation (euler-lotka formula)
% sum a (f(a)*I(a)*lambda^(-a-1)) = 1, where I(a) = p(a)*p(a-1)*...p(0)
% n(0,t) = n(0,0)*lambda^(t)

% data here
Amax = 2;
tmax = 100;
p = [0.5,0.25,0]; I = [1,0.5,0.125];
f = [0,1,5];
n = NaN(Amax,tmax);
n(1,1) = 100; n(2,1) = 50; n(3,1) = 13;

lambda_left = 0.25 ; % Euler-Lotka function -> +∞ as lambda -> 0
lambda_right = 5 ;   % Euler-Lotka function -> -1 as lambda -> +∞
lambdas = linspace(lambda_left,lambda_right,1000) ; % range of lambda
values = zeros(1,length(lambdas));

% plot to show
for j=1:length(lambdas)
    values(j) = euler_lotka_formula(lambdas(j),I,f,Amax); 
end
figure;
plot(lambdas,values,'LineWidth',1);
hold on;
plot(lambdas,zeros(size(values)),'r','LineWidth',1);

% compute zero crossing using matlab fzero
lambdabar = fzero(@(lambdas) euler_lotka_formula(lambdas,I,f,Amax),1);

n(1,:) = n(1,1)*lambdabar.^(1:tmax);
for a = 2:Amax+1
    n(a,:) = n(a-1,:).*p(a-1);
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

% euler-lotka formula
function value = euler_lotka_formula(lambda,I,f,Amax)
value = 0;
age = 1:Amax+1;
value = value + sum(f.*I.*lambda.^-age);
value = value - 1;
end