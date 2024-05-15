function mtg_check(rates, initial_cond, sigma, T, N,dt,flag)
% Martingality check
% INPUTS
% rates: interest rates
% initial_cond: initial value of the equity or property
% sigma: volatility
% T: time to maturity
% N: number of time steps
% dt: time grid
%
% OUTPUT
% plot of the Martingality check

% equity or property simulation
Simulated = simulate_GBM(rates, initial_cond, sigma, T, N); % simulations with no regular deduction
initial = initial_cond * ones(T)'; % vector of the initial value
Mean_sim = mean(Simulated(:,2:end),1)'; % vector of means at each time step
mean_discounted = exp(-rates.*dt).*Mean_sim;
deterministic = initial_cond*exp(rates.*dt)'; % vector of the deterministic projection

% plot of the Martinality check
figure
plot(dt, deterministic,'LineWidth',2)
hold on
plot(dt, Mean_sim,'LineWidth',2)
plot(dt, initial,'o')
plot(dt,mean_discounted,'*')
legend('Deterministic projection','Mean of stoch sims ','PF_{initial}','Discounted mean of stoch sims')

if flag == 1
    title('equity Martingality check')
    legend('Deterministic projection','Mean of stoch sims ','S_{initial}','Discounted mean of stoch sims')
else 
    title('property Martingality check')
    legend('Deterministic projection','Mean of stoch sims ','PF_{initial}','Discounted mean of stoch sims')
end
hold off

end % function mtg_check

