%% Logistic map cobwebbing

% for r < 1 try xinit = 0.2 0.5 0.8
% for 1<r<2 try xinit = 0.2 0.5 0.8
% for 2<r<3 try xinit = 0.3 0.7
r = 2.5;
f = @(x) r.*x.*(1-x);

figure
x_temp = 0 :0.001:1;
h1 = plot(x_temp,f(x_temp),'k','LineWidth',1);
hold on
h2 = plot(x_temp,x_temp,'k','LineWidth',1);

xinit = 0.78;
tmax = 100;
x = NaN(1,tmax);x(1) = xinit;
hold on
plot([xinit,xinit],[0,f(xinit)],'r','LineWidth',1);
for t = 2:tmax
    x(t) = f(x(t-1));
    hold on;
    plot([x(t-1),x(t)],[x(t),x(t)],'r','LineWidth',1);
    hold on
    plot([x(t),x(t)],[x(t),f(x(t))],'r','LineWidth',1);

end