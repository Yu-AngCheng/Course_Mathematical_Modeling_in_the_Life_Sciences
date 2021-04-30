%% EG Exercise 3.6.

clear
tmax = 1000000;
states = NaN(1,tmax);  

% Left stochastic matrix where aij means state j transfers to state i 
% Each column sum up to 1
A = [0.98, 0.10, 0.00;
     0.02, 0.70, 0.05;
     0.00, 0.20, 0.95;];

%% (a)
% Initial state
rng default;
Intial_state = 1;
states(1) = Intial_state;
random_number = rand(1,tmax);
for k= 1:tmax-1
    if random_number(k) < A(1,states(k))
        states(k+1) = 1;
    elseif random_number(k) < A(1,states(k))+A(2,states(k))
        states(k+1) = 2;
    else
        states(k+1) = 3;
    end
end
%% (b)
[V,D] = eig(A);D =diag(D);
index=find(abs(D-1)<1e-3); % to avoid the numerical flaws
dominant_eigvector = V(:,index);

% compare the observed occupancy with the equilibrium probability
state1_probability = sum(states==1)./tmax;
state2_probability = sum(states==2)./tmax;
state3_probability = sum(states==3)./tmax;
equilibrium_probability = dominant_eigvector./sum(dominant_eigvector);

%% (c)
% convert the first two indistinguishable states into one 
rstates = NaN(size(states));
rstates(states <= 2) = 1; % closed states
rstates(states == 3) = 2; % open states

%% (d)(e)
dwell_time_state1 = [];
dwell_time_state2 = [];
rstates = [rstates,999];
time_temp = 1;
for i = 2:tmax+1
    if rstates(i) == rstates(i-1)
        time_temp = 1 + time_temp;
    else
        if rstates(i-1) == 1
            dwell_time_state1 = [dwell_time_state1,time_temp];
        else
            dwell_time_state2 = [dwell_time_state2,time_temp];
        end
        time_temp = 1;
    end
end

nRun = 20;
numbins= 60;
figure
subplot(1,2,1);
h1 = histogram(dwell_time_state1,numbins,'Normalization','probability');
counts1 = h1.Values;
xvals1 = (h1.BinEdges(2:end)+h1.BinEdges(1:end-1))/2;
% one exponential
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf],...
                    'Upper', [Inf, 0],...
                    'StartPoint', [rand(), -rand()]);
    ft = fittype('a*exp(b*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals1', counts1', ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_close_one = fitObj_tmp{idx_opt};
a1_close_est_one = fitObj_close_one.a; lambda1_close_est_one = exp(fitObj_close_one.b);
% two exponential
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf, 0, -Inf],...
                    'Upper', [Inf, 0, Inf, 0],...
                    'StartPoint', [rand(), -rand(),rand(), -rand()]);
    ft = fittype('a*exp(b*x)+c*exp(d*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals1', counts1', ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_close_two = fitObj_tmp{idx_opt};
a1_close_est_two = fitObj_close_two.a; lambda1_close_est_two = exp(fitObj_close_two.b);
a2_close_est_two = fitObj_close_two.c; lambda2_close_est_two = exp(fitObj_close_two.d);
% three exponential(to exclude the effect of parameter numbers)
% the check the value of a3
% if a3 is near 0, then the third exponential is not necessary
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf, 0, -Inf, 0, -Inf],...
                    'Upper', [Inf, 0, Inf, 0, Inf, 0],...
                    'StartPoint', [rand(), -rand(), rand(), -rand(), rand(), -rand()]);
    ft = fittype('a*exp(b*x)+c*exp(d*x)+e*exp(f*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals1', counts1', ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_close_three = fitObj_tmp{idx_opt};
a1_close_est_three = fitObj_close_three.a; lambda1_close_est_three = exp(fitObj_close_three.b);
a2_close_est_three = fitObj_close_three.c; lambda2_close_est_three = exp(fitObj_close_three.d);
a3_close_est_three = fitObj_close_three.e; lambda3_close_est_three = exp(fitObj_close_three.f);
hold on; line1 = plot(fitObj_close_one,'-r',xvals1,counts1,'-g');
hold on; line1 = [line1, plot(fitObj_close_two,'-k',xvals1,counts1,'-g')];
hold on; line1 = [line1, plot(fitObj_close_three,'-b',xvals1,counts1,'-g')];
legend(line1([1 2 4 6]),'observed distribution curve','one exponential curve',...
   'two exponential curves','three exponential curves');
xlabel('Dwell times of closed states','Fontsize',12); ylabel('Normlized frequency','Fontsize',12);
legend boxoff; box off;
subplot(1,2,2);
h2 = histogram(dwell_time_state2,numbins,'Normalization','probability');
counts2 = h2.Values;
xvals2 = (h2.BinEdges(2:end)+h2.BinEdges(1:end-1))/2;
% one exponential
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf],...
                    'Upper', [Inf, 0],...
                    'StartPoint', [rand(), -rand()]);
    ft = fittype('a*exp(b*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals2', counts2',ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_open_one = fitObj_tmp{idx_opt};
a1_open_est_one = fitObj_open_one.a; lambda1_open_est_one = exp(fitObj_open_one.b);
% two exponential
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf, 0, -Inf],...
                    'Upper', [Inf, 0, Inf, 0],...
                    'StartPoint', [rand(), -rand(),rand(), -rand()]);
    ft = fittype('a*exp(b*x)+c*exp(d*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals2', counts2', ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_open_two = fitObj_tmp{idx_opt};
a1_open_est_two = fitObj_open_two.a; lambda1_open_est_two = exp(fitObj_open_two.b);
a2_open_est_two = fitObj_open_two.c; lambda2_open_est_two = exp(fitObj_open_two.d);
% three exponential(to exclude the effect of parameter numbers)
% the check the value of a3
% if a3 is near 0, then the third exponential is not necessary
for iRun = 1:nRun
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Lower', [0, -Inf, 0, -Inf, 0, -Inf],...
                    'Upper', [Inf, 0, Inf, 0, Inf, 0],...
                    'StartPoint', [rand(), -rand(), rand(), -rand(), rand(), -rand()]);
    ft = fittype('a*exp(b*x)+c*exp(d*x)+e*exp(f*x)', 'options', fo, 'independent', 'x');
    [fitObj_tmp{iRun}, gof] = fit(xvals2', counts2', ft);
    rmse_temp(iRun) = gof.rmse;
end
[~, idx_opt] = min(rmse_temp);fitObj_open_three = fitObj_tmp{idx_opt};
a1_open_est_three = fitObj_open_three.a; lambda1_open_est_three = exp(fitObj_open_three.b);
a2_open_est_three = fitObj_open_three.c; lambda2_open_est_three = exp(fitObj_open_three.d);
a3_open_est_three = fitObj_open_three.e; lambda3_open_est_three = exp(fitObj_open_three.f);
hold on; line2 = plot(fitObj_open_one,'r',xvals2,counts2,'-g');
hold on; line2 = [line2, plot(fitObj_open_two,'k',xvals2,counts2,'-g')];
hold on; line2 = [line2, plot(fitObj_open_three,'b',xvals2,counts2,'-g')];
legend(line2([1 2 4 6]),'observed distribution curve','one exponential curve',...
    'two exponential curves','three exponential curves')
xlabel('Dwell times of open states','Fontsize',12); ylabel('Normlized frequency','Fontsize',12);
legend boxoff; box off;