function [liabilities_rates , Lapse_BEL_rates, Death_BEL_rates, Expenses_BEL_rates, Commissions_BEL_rates, BOF_rates,delta_BOF_rates] = stress_rates(F0, ...
    P_death, lt, regular_deduction, COMM, expenses,dt, F_rates, benefit_commission, T,rates,BOF)


% calculate discouts
discounts = exp(-rates.*dt);

% computation of the Liabilities
[liabilities_rates , Lapse_BEL_rates, Death_BEL_rates, Expenses_BEL_rates, Commissions_BEL_rates] = Liabilities(F0, ...
          P_death, lt, regular_deduction, COMM, discounts, expenses,dt, F_rates, benefit_commission, T);

% computation of the BOF and deltaBOF with rates up
BOF_rates = F0 - liabilities_rates;
delta_BOF_rates = max(BOF-BOF_rates,0);

end