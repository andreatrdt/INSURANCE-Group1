function mtg_check(rates, S0, sigma_equity, T, N,dt)

% equity
S_equity = simulate_GBM(rates, S0, sigma_equity, T, N); % simulations with no regular deduction
S_initial = S0 * ones(T)'; % vector of the initial value
S_mean = mean(S_equity(:,2:end),1)'; % vector of means at each time step
S_mean_discounted = exp(-rates.*dt).*S_mean;
S_deterministic = S0*exp(rates.*dt)'; % vector of the deterministic projection

% plot of the Martinality check for Equity
figure
plot(dt, S_deterministic,'LineWidth',2)
hold on
plot(dt, S_mean,'LineWidth',2)
plot(dt, S_initial,'o')
plot(dt,S_mean_discounted,'*')
legend('Deterministic projection','Mean of stoch sims ','PF_{initial}','Discounted mean of stoch sims')
title('Equity Martingality check')
hold off

end