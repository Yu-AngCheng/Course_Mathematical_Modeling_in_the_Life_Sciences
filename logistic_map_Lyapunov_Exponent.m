%% logistic map Lyapunov Exponent
clear

dr=0.001;
rmin = 0;
rmax = 4;
rr=rmin:dr:rmax;
tmax = 10000;
x0 = 0.6;
Lyapunov_Exponent = NaN(1,length(rr));
figure
for i = 1:length(rr)
    r = rr(i);
    f = @(x) r.*x.*(1-x); % the function f
    x_temp = NaN(1,tmax);
    x_temp(1) = x0;
    for t = 2:tmax
        x_temp(t) = f(x_temp(t-1));
    end
    Lyapunov_Exponent(i) = 1/tmax*sum(log(abs(r-2*r*x_temp)));
end
plot(rr,Lyapunov_Exponent);
ylabel('Lyapunov Exponent','FontSize',15);xlabel('r','FontSize',15)