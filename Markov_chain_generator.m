%% Markov_chain_generator

tmax = 100;
statet = NaN(1,tmax);  

% Left stochastic matrix where aij means state j transfers to state i 
% Each column sum up to 1
A = [0.75, 0.15 ;
     0.25, 0.85] ;
% Initial state
statet(1) = 1;

random_number = rand(1,tmax);
for k= 1:tmax-1
    if random_number(k) < A(1,statet(k))
        statet(k+1) = 1;
    else
        statet(k+1) = 2;
    end
end

figure;
plot(statet,'o-')
xlabel('timestep','FontSize',16)
ylabel('state','FontSize',16)