%% EG Lab Manual 11.2, including the Matlab challenge
% Convert the sequence of 1000 random numbers r from the previous exercise
% into a sequence of outcomes of coin tosses in which the probability of 
% heads is 0:6 and the probability of tails is 0:4. 
% Let 1 represent an outcome of heads and let 0 represent an outcome
% of tails.

clear

% Generate a sequence r of 1000 random numbers uniformly distributed 
% in the unit interval [0; 1]

num = 1000;
rng default; r = rand(1,num);

% using for loop
seq = zeros(1,num);
for i = 1:num
    if r(i) < 0
        seq(i) = 1;
    end
end

% without using for loop
seq = NaN(1,num);
seq(r < 0.6) = 1; %  head
seq(r >= 0.6) = 0; % tail

% Show the binomial distribution in the coin tossing experiemnt
% when n is large, the binomial distribution converges to normal
% distribution

% the analytic solution
head_probability = 0.6;
tail_probability = 1  - head_probability;
num_trial = 1000;
k_interest = 500:1:700;
for i = 1:length(k_interest)
    k = k_interest(i);
    probability_k(i) = nchoosek(num_trial,k).*head_probability.^k.*...
    tail_probability.^(num_trial - k);
end
figure
plot(k_interest,probability_k,'-k','LineWidth',1);
xlabel('# of heads','FontSize',12);
ylabel('probability','FontSize',12);
box off

% the experiment 1
num_experiement1 = 1000; num_trial1 = 1000;
rng default;
rr1 = rand(num_experiement1,num_trial1);

seq1 = NaN(size(rr1));
seq1(rr1 < head_probability) = 1; %  head
seq1(rr1 >= head_probability) = 0; % tail

Number_of_heads1 = sum(seq1,2);
figure;
xlabel('# of heads','FontSize',12);
yyaxis left;
histogram(Number_of_heads1); 
ylabel('Frequency','FontSize',12);
hold on;
yyaxis right
plot(k_interest,probability_k,'-k','LineWidth',1); 
ylabel('Probability','FontSize',12);

% the experiment 2
num_experiement2 = 10000; num_trial2 = 100;
rng default;
rr2 = rand(num_experiement2,num_trial2);

seq2 = NaN(size(rr2));
seq2(rr2 < head_probability) = 1; %  head
seq2(rr2 >= head_probability) = 0; % tail

Number_of_heads2 = sum(seq2,2);

% to compare the results between 2 experiment, we zscore the both and draw
% it in the normlized version.
% Compare the two distributions with the standard normal distribution
figure
binwidth = 0.12;
edge = -4:binwidth:4;
histogram(zscore(Number_of_heads1),edge,'Normalization','probability');
hold on;
histogram(zscore(Number_of_heads2),edge,'Normalization','probability');
hold on
plot(edge,normpdf(edge)*binwidth,'-k','LineWidth',1.5);
legend({'Experiment 1','Experiment 2',...
    'Normal Distribution'})
legend boxoff; box off;
xlabel('zscore of # of heads','FontSize',12);
ylabel('Normalized Frequency / Probability')
