function [S] = simulate_GBM(S0, sigma, T, N, M , regular_deduction)
    % INPUT
    % S0: initial stock price
    % mu: expected return
    % sigma: volatility
    % T: time to maturity
    % N: number of time steps
    % M: number of paths
    % OUTPUT
    % S: stock price paths
    % monte carlo simulation

    dt = [1; T(2:end, 1)-T(1:end-1, 1)]''; 

    S = zeros(N, length(dt)); 
    S(:,1) = S0*ones(N,1);

    for i = 1:(length(dt)-1)
    S(:,i+1) = (1-regular_deduction)*S(:,i).*exp((-(sigma^2)/2)*dt(i)+sigma*randn(N,1)); 
    end

    
end