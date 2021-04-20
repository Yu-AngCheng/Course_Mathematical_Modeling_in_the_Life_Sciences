%% Logistic map Bifurcation
clear

dr=0.001;
rmin = 0;
rmax = 4;
rr=rmin:dr:rmax;
tmax = 1000;
x0 = 0.5;
ncmax = 100; % choose the last ncmax iterations and assume the convergence
figure
for r = rr
    f = @(x) r.*x.*(1-x); % the function f
    x_temp = NaN(1,tmax);
    x_temp(1) = x0;
    for t = 2:tmax
        x_temp(t) = f(x_temp(t-1));
    end
    hold on
    plot(r*ones(1,ncmax),x_temp(1,end-ncmax+1:end),'ro','MarkerSize',0.5);
end
ylabel('attractor','FontSize',15);xlabel('r','FontSize',15)