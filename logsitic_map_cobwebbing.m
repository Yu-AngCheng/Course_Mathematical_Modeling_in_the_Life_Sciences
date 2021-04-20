%% Logistic map cobwebbing
% What we can analyze is what we can sketch!
% Just sketch the process of iterations 

% Short summary:
% for 0 < r < 1, we have only one stable steady solution xbar = 0
% and f'(xbar) = r ∈(0,1)
% for 1 < r < 3, we also have only one stable steady solution xbar = 1-1/r
% and f'(xbar) = 2 - r. So if r > 2, f'(xbar) < 0, causing spirals
% when r < 3 it spirals in, but when r > 3, it spirals out.


% for 0 < r < 1 try xinit = 0.2 0.5 0.8
% for 1 < r < 2 try xinit = 0.2 0.5 0.8
% for 2 < r < 3 try xinit = 0.3 0.7
% for r > 3, just try r = 3.1, 3.45, 3.56

clear

r = 1;
f = @(x) r.*x.*(1-x);

figure
x_temp = 0 :0.001:1;
h1 = plot(x_temp,f(x_temp),'k','LineWidth',1); % the function f
hold on
h2 = plot(x_temp,x_temp,'k','LineWidth',1); % the diagonal line

xinit = 0.9;
tmax = 100;
x = NaN(1,tmax);x(1) = xinit;
hold on
plot([xinit,xinit],[0,f(xinit)],'Color','b','LineWidth',0.5);
for t = 2:tmax
    x(t) = f(x(t-1));
    if (tmax-t)<40 % the last 40 iterations
        hold on;
        plot([x(t-1),x(t)],[x(t),x(t)],'-ro','LineWidth',1.5);
        hold on
        plot([x(t),x(t)],[x(t),f(x(t))],'-ro','LineWidth',1.5);
    else
        hold on;
        plot([x(t-1),x(t)],[x(t),x(t)],'b','LineWidth',0.5);
        hold on
        plot([x(t),x(t)],[x(t),f(x(t))],'b','LineWidth',0.5);
    end
end
axis square

% show the process
figure
plot(x);