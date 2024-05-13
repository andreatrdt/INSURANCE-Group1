function [liabilities_shocked_pr, Lapse_BEL_property, Death_BEL_property, Expenses_BEL_property, Commissions_BEL_property,delta_BOF_pr] = ...
     property_risk(S0, PF_0, rates, sigma_pf, T, N, regular_deduction, P_death, lt, COMM, discounts, expenses, dt, benefit_commission, BOF,S,F0)

% Initial value of the shocked value of the property
P0_shocked=(1-0.25)*PF_0;

% initial value for the fund
F0_property = P0_shocked + S0;  

% property simulation
P_shocked = simulate_GBM(rates(1:T), P0_shocked, sigma_pf, T, N, regular_deduction);

% value of the fund 
F = P_shocked + S;   

% computation of Liabilities 
[liabilities_shocked_pr, Lapse_BEL_property, Death_BEL_property, Expenses_BEL_property, Commissions_BEL_property] = Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts, expenses, dt, F, benefit_commission, T);

% computation of BOF and delta BOF in case of property risk
BOF_pr = F0_property - liabilities_shocked_pr;
delta_BOF_pr = max(BOF - BOF_pr,0);

end