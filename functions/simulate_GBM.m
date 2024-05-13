function S = simulate_GBM(rates, S0, sigma, T, N, regular_deduction)
    % Simulates GBM
    % INPUTS
    % S0: initial stock price
    % mu: expected return
    % sigma: volatility
    % T: time to maturity
    % N: number of time steps
    %
    % OUTPUT
    % S: simulation of a Geometric Brownian Motion
    
    % if we do not give regular_deduction as an input, then we are in the
    % Martingality Check and RD has to be set to 0

    % fix random seed 
    Var_seed = 42; % the answer to the ultimate question of life, the universe, and everything
    rng(Var_seed)

    if nargin <= 5
        regular_deduction = 0;
    end

    % time grid
    dt = (1:1:T)';
    % time step
    delta_time = dt - [0; dt(1:end-1)];

    % initialization of the variables
    S = zeros(N,T); 
    S(:,1) = S0*ones(N,1);

    % computation of the discounts
    discounts = exp(-rates.*dt);

    % computation of the forward rates
    fwd = discounts./[1; discounts(1:end-1)];
    fwd_rates = -log(fwd);

    % random variable ~ N(0,1)
    g = randn(N,length(dt));
    
    for i = 1:length(dt)
        S(:,i+1) = (1-regular_deduction)*S(:,i).*exp((fwd_rates(i)-(sigma^2)/2)*delta_time(i)+sigma*sqrt(delta_time(i))*g(:,i)); 
    end
    
end % function simulate_GBM







