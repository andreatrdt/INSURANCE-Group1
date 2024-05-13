function [liabilities_shocked_eq, Lapse_BEL_equity, Death_BEL_equity, Expenses_BEL_equity, Commissions_BEL_equity,delta_BOF_eq] = ...
    equity_risk(S0, PF_0, rates, sigma_equity, T, N, regular_deduction, P_death, lt, COMM, discounts, expenses, dt, PF, benefit_commission, BOF,F0)

symm_adj = 0.0525; % march 2024 from EIOPA

% initial value for shocked equity
S0_shocked = (1 - 0.39 - symm_adj) * S0; % 0.39 since it is of type 1

% initial value for the fund
F0_equity = S0_shocked + PF_0;

% equity simulation
S_shocked = simulate_GBM(rates(1:T), S0_shocked, sigma_equity, T, N, regular_deduction);

% value of the fund
F = S_shocked + PF; 

[liabilities_shocked_eq, Lapse_BEL_equity, Death_BEL_equity, Expenses_BEL_equity, Commissions_BEL_equity] = Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

BOF_eq = F0_equity - liabilities_shocked_eq;
delta_BOF_eq = max(BOF - BOF_eq,0);

end