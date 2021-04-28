%% EG Lab Manual 11.1
% Coin tossing Exercise 11.1
% Experiment with sequences of coin flips produced by a random number 
% generator

clear

% Generate a sequence r of 1000 random numbers uniformly distributed 
% in the unit interval [0; 1]

num = 1000;
rng default; r = rand(1,num);

% draw the histogram and compare the observed value to the expected value
nbins = 10; edge = 0:0.1:1;
figure
h = histogram(r,edge);
hold on; yline(num/nbins,'k--','LineWidth',1);
xlabel('r','FontSize',12);ylabel('Counts','FontSize',12)
legend({'observed value','expected value'},'Orientation',...
    'horizontal','Location','SouthOutside');
legend boxoff; box off;

% to further verify the uniform distribution, do descriptive statistics and
% hypothesis testing
sample_mean = mean(r); sample_var = var(r);
counts = h.Values;
chi2 = sum((counts-num/nbins).^2./(num/nbins));
p = chi2cdf(chi2,nbins-1,'upper');


% --------------------------------------------------------------------%
% to check the CLT theorm
seq_length = 100000;
seq_num = 1000;
rng default 
rr = rand(seq_num,seq_length); rr_mean = mean(rr,2);
% the sample mean follows a normal distribution with mean = 0.5 and
% variance = (1/12)/10000;

% to do hypothesis testing to test the CLT theorm
analytic_mean = 0.5;
analytic_variance = (1/12)/seq_length;
[h,p] = kstest((rr_mean - analytic_mean)./sqrt(analytic_variance));
E_rr_mean = mean(rr_mean);
E_rr_var = var(rr_mean);


% --------------------------------------------------------------------%
% to check the CLT theorm
seq_length = 10000;
seq_num = 1000;
rng default 
rr = rand(seq_num,seq_length); rr_mean = mean(rr,2);
% the sample mean follows a normal distribution with mean = 0.5 and
% variance = (1/12)/10000;

% to do hypothesis testing to test the CLT theorm
analytic_mean = 0.5;
analytic_variance = (1/12)/seq_length;
[h,p] = kstest((rr_mean - analytic_mean)./sqrt(analytic_variance));
E_rr_mean = mean(rr_mean);
E_rr_var = var(rr_mean);