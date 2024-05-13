function [liabilities_shocked_mortality, Lapse_BEL_mortality, Death_BEL_mortality, Expenses_BEL_mortality, Commissions_BEL_mortality,delta_BOF_mortality] = ...
    mortality_risk(F0, T, regular_deduction, P_death, lt, COMM, discounts, expenses, dt, benefit_commission, BOF , S ,PF)

P_death_shocked = min(1, P_death*1.15);

% simulate equity prices
F = S + PF;

% Computation of Liabilities
[liabilities_shocked_mortality, Lapse_BEL_mortality, Death_BEL_mortality, Expenses_BEL_mortality, Commissions_BEL_mortality] = Liabilities(F0, ...
            P_death_shocked, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% computation of BOF and delta BOF in case of mortality risk
BOF_mortality = F0 - liabilities_shocked_mortality;
delta_BOF_mortality = max(BOF - BOF_mortality,0);

end