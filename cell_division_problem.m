%% cell devision problem
% synchronously dividing cells
% each cell divides into 4
% Fraction p of cell survive

islinear = 0;

n0 = 100; % specify inital population
tmax = 1000; % specify final time
n = NaN(1,tmax); % Number of populations
n(1) = n0;

if islinear    
    % parameters that can be tried
    % p = 0.24;
    % p = 0.26;
%     p = 0.25;
    for t = 2:tmax
        n(t) = n(t-1) * 4 * p;
    end
else 
    % logistic map here
    % cannot be solved analytically
    % use numerical simulation!
    
    % parameters that can be tried
%     p = 0.24;
%     p = 0.26;
    p = 0.3;
%     p = 0.75;
%     p = 0.90;
    K = 1000;
    for t = 2:tmax
        n(t) = n(t-1) * 4 * p * (1 - n(t-1)/K);
    end
end

figure
plot(1:tmax,n);