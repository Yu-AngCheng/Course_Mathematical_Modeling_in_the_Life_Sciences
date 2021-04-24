%% EG exercise 3.5


clear


tmax = 100000;
statet = NaN(1,tmax);  

% Left stochastic matrix where aij means state j transfers to state i 
% Each column sum up to 1
A = [0.75, 0.15 ;
     0.25, 0.85] ;
dwell_time = 1:1:50;
dwell_time_1_analytic = A(1,1).^(dwell_time-1).*(1-A(1,1));
dwell_time_2_analytic = A(2,2).^(dwell_time-1).*(1-A(2,2));

Intial_state = 1;
statet(1) = Intial_state;
rng default;
random_number = rand(1,tmax);
for k= 1:tmax-1
    if random_number(k) < A(1,statet(k))
        statet(k+1) = 1;
    else
        statet(k+1) = 2;
    end
end

dwell_time_state1 = [];
dwell_time_state2 = [];
% add one to avoid the last
statet = [statet,999];
time_temp = 1;
for i = 2:tmax+1
    if statet(i) == statet(i-1)
        time_temp = 1 + time_temp;
    else
        if statet(i-1) == 1
            dwell_time_state1 = [dwell_time_state1,time_temp];
        else
            dwell_time_state2 = [dwell_time_state2,time_temp];
        end
        time_temp = 1;
    end
end

% compare the observed distributions of dwell times and the analytic ones.
figure
subplot(1,2,1);
histogram(dwell_time_state1,'Normalization','probability');
xlabel('Dwell times','Fontsize',12); ylabel('Normlized frequency','Fontsize',12);
hold on; plot(dwell_time,dwell_time_1_analytic,'--','LineWidth',1.5);
legend({'observed distribution','analytic distribution'})
legend boxoff; box off;
subplot(1,2,2);
histogram(dwell_time_state2,'Normalization','probability');
xlabel('Dwell times','Fontsize',12); ylabel('Normlized frequency','Fontsize',12);
hold on; plot(dwell_time,dwell_time_2_analytic,'--','LineWidth',1.5);
legend({'observed distribution','analytic distribution'})
legend boxoff; box off;
