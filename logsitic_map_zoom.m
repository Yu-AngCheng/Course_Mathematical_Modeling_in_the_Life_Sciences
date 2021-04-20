%% Logistic map zoom
% plot f(x) vs. x, f(f(x)) vs. x and f(f(f(f(x))))

xx=0:0.001:1;
r=(1+sqrt(6)+0.01);
f = @(x) r.*x.*(1-x); % the function f

fx=f(xx);
ffx=f(fx);

figure;
subplot(1,2,1);
plot(xx,fx,'LineWidth',2);
xlabel('x','FontSize',15);
ylabel('f(x)','FontSize',15);
axis equal;
axis([0 1 0 1]);
hold on;
plot(xx,xx,'r--')

subplot(1,2,2);
plot(xx,ffx,'LineWidth',2);
xlabel('x','FontSize',15);
ylabel('f(f(x))','FontSize',15);
axis equal;
axis([0 1 0 1]);
hold on;
plot(xx,xx,'r--')

