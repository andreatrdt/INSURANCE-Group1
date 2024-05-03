function [S] = simulate_GBM(rates, S0, sigma, T, N, regular_deduction)
    % INPUT
    % S0: initial stock price
    % mu: expected return
    % sigma: volatility
    % T: time to maturity
    % N: number of time steps
    % OUTPUT
    % S: stock price paths
    
    dt = [1:1:T]'; 
    S = zeros(N,T); 
    S(:,1) = S0*ones(N,1);
    %calculate discouts
    discounts = exp(-rates.*dt);
    %calculate forward rates
    fwd = discounts./[1; discounts(1:end-1)];
    fwd_rates = -log(fwd);
    delta_time = dt - [0; dt(1:end-1)];
    g = randn(N,length(dt));
    
    for i = 1:length(dt)
    S(:,i+1) = (1-regular_deduction)*S(:,i).*exp((fwd_rates(i)-(sigma^2)/2)*delta_time(i)+sigma*sqrt(delta_time(i))*g(:,i)); 
    end
    
    
end