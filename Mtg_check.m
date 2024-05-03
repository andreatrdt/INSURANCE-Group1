function check = Mtg_check(S,rates,dt,S0,discounts)
% INPUT
% S: stock price paths
% rates: interest rates
% dt: time steps
% S0: initial stock price
% OUTPUT
% check: plot of the stock price paths



S_mean = mean(S,1)';
dt= [0 ; dt ];
rates = [0, rates];

figure;
plot(dt,S_mean)
hold on
plot(dt,S0*exp(rates.*dt),'r')

check = 1;

end
