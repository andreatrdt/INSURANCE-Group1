function [check] = Mtg_check(S,dt,S0,RD)
% INPUT
% S: stock price paths
% rates: interest rates
% dt: time steps
% S0: initial stock price
% OUTPUT
% check: plot of the stock price paths

% check if the stock price paths are martingales
for i = 1:length(dt)
    S_mean(i) = mean(S(:,i));
end


figure 
plot(dt,S_mean)
hold on
plot(dt,(1-RD)*S0*ones(length(dt)),'r')
hold off
xlabel('Time')
legend('Simulation ','Expected simulation','50 years expected simulation')

check = 1;
end
