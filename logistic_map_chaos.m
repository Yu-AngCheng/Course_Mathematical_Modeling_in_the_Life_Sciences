clear

r = 4;
f = @(x) r.*x.*(1-x);
tmax = 1000000;
figure
x_temp = 0 :0.001:1;
h1 = plot(x_temp,f(x_temp),'k','LineWidth',1); % the function f
hold on
h2 = plot(x_temp,x_temp,'k','LineWidth',1); % the diagonal line

xinit = 0.6;
x1 = NaN(1,tmax);x1(1) = xinit;
hold on
plot([xinit,xinit],[0,f(xinit)],'Color','b','LineWidth',0.5);
for t = 2:tmax
    x1(t) = f(x1(t-1));
    hold on;
    plot([x1(t-1),x1(t)],[x1(t),x1(t)],'b','LineWidth',0.5);
    hold on
    plot([x1(t),x1(t)],[x1(t),f(x1(t))],'b','LineWidth',0.5);
end

xinit = 0.6000001;
x2 = NaN(1,tmax);x2(1) = xinit;
hold on
plot([xinit,xinit],[0,f(xinit)],'Color','r','LineWidth',0.5);
for t = 2:tmax
    x2(t) = f(x2(t-1));
    hold on;
    plot([x2(t-1),x2(t)],[x2(t),x2(t)],'r','LineWidth',0.5);
    hold on
    plot([x2(t),x2(t)],[x2(t),f(x2(t))],'r','LineWidth',0.5);
end
axis square

% show the process
% x1 and x2 only share the same points only at early iterations and differ
% at later iterations
figure
plot(x1);
hold on
plot(x2);
xlabel('Iterations');
ylabel('x_n');

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
