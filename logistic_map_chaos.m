%% show the chaos 
clear

r = 4;
f = @(x) r.*x.*(1-x);
tmax = 10000;
xinit = 0.6;
x1 = NaN(1,tmax);x1(1) = xinit;
xinit = 0.6000001;
x2 = NaN(1,tmax);x2(1) = xinit;
for t = 2:tmax
    x1(t) = f(x1(t-1));
    x2(t) = f(x2(t-1));
end

% show the process
% x1 and x2 only share the same points only at early iterations and differ
% at later iterations
figure
plot(x1);
hold on
plot(x2);
xlabel('Iterations');
ylabel('x_n');
xlim([0,50])

% show the distribution
% But x1 and x2 have the same distribution
% For better illustration, choose tmax =10000;
figure
subplot(3,1,1);
histogram(x1);
subplot(3,1,2);
histogram(x2);
subplot(3,1,3);
histogram(x1);
hold on;
histogram(x2);