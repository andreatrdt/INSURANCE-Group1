%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment Insurance
% Authors:
%   - Giovanni Riondato 
%   - Giacomo Manfredi
%   - Lorenzo Tolomelli
%   - Andera Tarditi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial settings
clc;
clear all;
close all;

Var_seed = 42; % the answer to the ultimate question of life, the universe, and everything
rng(Var_seed)

% select format
format bank
format long

% fix random seed

% start run time
tic

% load data from xlsx file

% filename in .xls
filename = 'Life_tables_of_the_resident_population.xlsx';
 
% read excel data from Life_tables_of_the_resident_population.xlsx
P_death = xlsread(filename,1,'D63:D112'); % since the client is a 60 years old male

P_death = P_death./1000;

% read excel data from EIOPA
filename = 'EIOPA_RFR_20240331_Term_Structures.xlsx';

rates_UP = xlsread(filename,5,'S11:S160');
rates_DOWN = xlsread(filename,6,'S11:S160');
rates = xlsread(filename,1,'S11:S160');

% convert discrete rates to continuous rates
rates_UP = log(1+rates_UP);
rates_DOWN = log(1+rates_DOWN);
rates = log(1+rates);

% plot yield rate curve
figure
hold on 
years = 1:150;
plot(years,rates)
plot(years,rates_UP)
plot(years,rates_DOWN)
rate_50 = interp1(years,rates,50);
plot(50,rate_50,'ro')
legend( 'rates','rates_{UP}','rates_{DOWN}')

%% Initialization of the parameters

% assets
F0 = 1e5; % value of fund at t = 0
S0 = 0.8*F0; % value of equity at t = 0
sigma_equity = 0.2; % volatility of equity
PF_0 = 0.2*F0; % value of property at t = 0
sigma_pf = 0.1; % volatility of property features

% liabilities

% contract terms
lapse_pay = 20; % benefit commission
regular_deduction = 0.022; % regular deduction
COMM = 0.014; % commission

% other specifications
T = 50; % number of years
dt = (1:1:T)';  % time grid

% operating & economic assumptions
lt = 0.15 * ones(length(dt),1); % lapse rate
inflation = 0.02; % inflation rate
yearly_expense_t0 = 50; % yearly expenses at t = 0
expenses = yearly_expense_t0.*(1+inflation).^[0; dt(1:end-1)]; % expenses

disp('Parameters loaded...')

% consider only the first T rates
rates = rates(1:T);
rates_UP = rates_UP(1:T);
rates_DOWN = rates_DOWN(1:T);

%% Martingality Check

% in order to choose a proper N, we want to check Martingality for Equity
% and Property

disp('Martingality check...')

N = 1e6; % we want to check whether this value is enough

mtg_check(rates, S0, sigma_equity, T, N, dt); % Equity
mtg_check(rates, PF_0, sigma_pf, T, N, dt); % Property

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intrest rate risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Interest rate risk analysis...')

%% Basic scenario

% simulate equity prices
S = simulate_GBM(rates, S0, sigma_equity, T, N, regular_deduction);

% simulate property features
PF = simulate_GBM(rates, PF_0, sigma_pf, T, N, regular_deduction);

% calculate fund value
F = S + PF;

% calculate discouts
discounts = exp(-rates.*dt);

% computation of the Liabilities
[liabilities, Lapse_BEL, Death_BEL, Expenses_BEL, Commissions_BEL] = Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts, expenses, dt, F, lapse_pay, T);

% computation of the BOF in basic scenario
BOF = F0 - liabilities;

% plot paths 
% figure
% hold on
% for i=1:N
%     plot(dt,PF_simulated(i,:))
% end
% figure
% hold on
% for i=1:N
%     plot(dt,S_simulated(i,:))
% end

%% Stress scenario UP

%simulate equity prices
S_rates_UP = simulate_GBM(rates_UP, S0, sigma_equity, T, N, regular_deduction);

% simulate property features
PF_rates_UP = simulate_GBM(rates_UP, PF_0, sigma_pf, T, N, regular_deduction);

% calculate fund value
F_rates_UP = S_rates_UP + PF_rates_UP;

% calculate discouts
discounts_UP = exp(-rates_UP.*dt);

% computation of the Liabilities
[liabilities_rates_UP , Lapse_BEL_rates_UP, Death_BEL_rates_UP, Expenses_BEL_rates_UP, Commissions_BEL_rates_UP] = Liabilities(F0, ...
          P_death, lt, regular_deduction, COMM, discounts_UP, expenses,dt, F_rates_UP, lapse_pay, T);

% computation of the BOF and deltaBOF with rates up
BOF_rates_UP = F0 - liabilities_rates_UP;
delta_BOF_rates_UP = max(BOF-BOF_rates_UP,0);

%% Stress scenario DOWN

% simulate equity prices
S_rates_DOWN = simulate_GBM(rates_DOWN, S0, sigma_equity, T, N, regular_deduction);

% simulate property features
PF_rates_DOWN = simulate_GBM(rates_DOWN, PF_0, sigma_pf, T, N, regular_deduction);

%calculate discouts
discounts_DOWN = exp(-rates_DOWN.*dt);

% calculate fund value
F_rates_DOWN = S_rates_DOWN + PF_rates_DOWN;

% computation of the Liabilities
[liabilities_rates_DOWN, Lapse_BEL_rates_DOWN, Death_BEL_rates_DOWN, Expenses_BEL_rates_DOWN, Commissions_BEL_rates_DOWN] = Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts_DOWN, expenses,dt,F_rates_DOWN, lapse_pay,T);

% computation of the BOF and deltaBOF with rates down
BOF_rates_DOWN = F0 - liabilities_rates_DOWN;
delta_BOF_rates_DOWN = max(BOF-BOF_rates_DOWN,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equity risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Equity risk analysis...')
% introduction of the symmetric adjustment
symm_adj = 0.0525; % march 2024 from EIOPA

% initial value for shocked equity
S0_shocked = (1 - 0.39 - symm_adj) * S0; % 0.39 since it is of type 1

% initial value for the fund
F0_equity = S0_shocked + PF_0;

% equity simulation
S_shocked = simulate_GBM(rates(1:T), S0_shocked, sigma_equity, T, N, regular_deduction);

% value of the fund
F = S_shocked + PF; 

% computation of Liabilities 
[liabilities_shocked_eq, Lapse_BEL_equity, Death_BEL_equity, Expenses_BEL_equity, Commissions_BEL_equity] = Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% computation of BOF and delta BOF in case of equity risk
BOF_eq = F0_equity - liabilities_shocked_eq;
delta_BOF_eq = max(BOF - BOF_eq,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Property risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Property risk analysis...')

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
            P_death, lt, regular_deduction, COMM, discounts, expenses, dt, F, lapse_pay, T);

% computation of BOF and delta BOF in case of property risk
BOF_pr = F0_property - liabilities_shocked_pr;
delta_BOF_pr = max(BOF - BOF_pr,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORTALITY RISK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Mortality risk analysis...')

P_death_shocked = min(1, P_death*1.15);

% simulate equity prices
F = S + PF;

% Computation of Liabilities
[liabilities_shocked_mortality, Lapse_BEL_mortality, Death_BEL_mortality, Expenses_BEL_mortality, Commissions_BEL_mortality] = Liabilities(F0, ...
            P_death_shocked, lt, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% computation of BOF and delta BOF in case of mortality risk
BOF_mortality = F0 - liabilities_shocked_mortality;
delta_BOF_mortality = max(BOF - BOF_mortality,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lapse risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Lapse risk analysis...')

% Lapse rate
lt = 0.15;

%% lapse UP risk

lt_shocked_UP = min(1.5*lt,1) * ones(length(dt),1);

% computation of the Liabilities
[liabilities_shocked_lapse_UP, Lapse_BEL_lapse_UP, Death_BEL_lapse_UP, Expenses_BEL_lapse_UP, Commissions_BEL_lapse_UP] = Liabilities(F0, ...
            P_death, lt_shocked_UP, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% computation of BOF and delta BOF
BOF_lapse_UP = F0 - liabilities_shocked_lapse_UP;
delta_BOF_lapse_UP = max(BOF - BOF_lapse_UP,0);

%% lapse DOWN risk

lt_shocked_DOWN = max(0.5*lt,lt-0.2) * ones(length(dt),1);

% Computation of the Liabilities
[liabilities_shocked_lapse_DOWN, Lapse_BEL_lapse_DOWN, Death_BEL_lapse_DOWN, Expenses_BEL_lapse_DOWN, Commissions_BEL_lapse_DOWN] = Liabilities(F0, ...
            P_death, lt_shocked_DOWN, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% computation of BOF and delta BOF
BOF_lapse_DOWN = F0 - liabilities_shocked_lapse_DOWN;
delta_BOF_lapse_DOWN = max(BOF - BOF_lapse_DOWN,0);



%% Lapse mass risk

lt_shocked_mass = [lt + 0.4;lt*ones(T-1,1)];     % verificare se + 0.4 0 direttamente 0.4

% Computation of the Liabilities
[liabilities_shocked_lapse_mass, Lapse_BEL_lapse_mass, Death_BEL_lapse_mass, Expenses_BEL_lapse_mass, Commissions_BEL_lapse_mass] = Liabilities(F0, ...
            P_death, lt_shocked_mass, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% computation of BOF and delta BOF
BOF_lapse_mass = F0 - liabilities_shocked_lapse_mass;
delta_BOF_lapse_mass = max(BOF - BOF_lapse_mass,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Expense risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Expense risk analysis...')

% reset lapse rate
lt = 0.15*ones(length(dt),1);

% compute the shocked expenses
expenses_t0_shocked = yearly_expense_t0*1.1;
expenses_shocked = expenses_t0_shocked.*(1+inflation+0.01).^[0; dt(1:end-1)];

% computation of the Liabilities
[liabilities_shocked_expense, Lapse_BEL_expense, Death_BEL_expense, Expenses_BEL_expense, Commissions_BEL_expense]= Liabilities(F0, ...
            P_death, lt, regular_deduction, COMM, discounts, expenses_shocked,dt,F,lapse_pay,T);

% computation of the BOF and delta BOF in case of expense risk
BOF_expense = F0 - liabilities_shocked_expense;
delta_BOF_expense = max(BOF - BOF_expense,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Catastrophe risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Catastrophe risk analysis...')

% Catastrophe scenario
P_death_cat = [P_death(1)+0.0015; P_death(2:end)];

% Computation of the Liabilities 
[liabilities_cat, Lapse_BEL_cat, Death_BEL_cat, Expenses_BEL_cat, Commissions_BEL_cat]= Liabilities(F0, ...
            P_death_cat, lt, regular_deduction, COMM, discounts, expenses,dt,F,lapse_pay,T);

% Delta BOF
BOF_cat = F0 - liabilities_cat;
delta_BOF_catastrophe = max(BOF - BOF_cat,0);

% Basic own fund and delta BOF
BOF_cat = F0 - liabilities_cat;
delta_BOF_cat = max(BOF-BOF_cat,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BSCR Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Basic Solvency Capital requirment computation...')

% SCR Market risk
SCR_IR = max(delta_BOF_rates_UP , delta_BOF_rates_DOWN); % Interest rate risk
SCR_EQ = delta_BOF_eq; % Equity risk
SCR_PR = delta_BOF_pr; % Property risk 

MKT_vec = [SCR_IR ; SCR_EQ ; SCR_PR]; 

% check if we are exposed to IR up or down to choose the correlation matrix
if delta_BOF_rates_UP > delta_BOF_rates_DOWN
    MKT_corr = [1 0 0; 0 1 0.75; 0 0.75 1];
else
    MKT_corr = [1 0.5 0.5; 0.5 1 0.75; 0.5 0.75 1];
end

SCR_MKT = sqrt(MKT_vec' * MKT_corr * MKT_vec);

% SCR Life risk
SCR_MORT = delta_BOF_mortality; % Mortality risk
SCR_LAPSE = max([delta_BOF_lapse_UP,delta_BOF_lapse_DOWN,delta_BOF_lapse_mass]); % Lapse risk
SCR_EXP = delta_BOF_expense; % Expense risk
SCR_CAT = delta_BOF_catastrophe; % Catastrophe risk
LIFE_vec = [SCR_MORT ; SCR_LAPSE ; SCR_EXP ; SCR_CAT];

% correlation matrix
LIFE_corr = [1 0 0.25 0.25; 0 1 0.5 0.25; 0.25 0.5 1 0.25 ; 0.25 0.25 0.25 1 ];

SCR_LIFE = sqrt(LIFE_vec' * LIFE_corr * LIFE_vec);

%% computation of BSCR

SCR = [SCR_MKT ; SCR_LIFE];

% correlation matrix
SCR_corr = [1 0.25; 0.25 1]; 

BSCR = sqrt(SCR' * SCR_corr * SCR);

% print results
fprintf('Results:\n');
fprintf('---------------------------------\n');
fprintf('BSCR:     %f\n', BSCR);
fprintf('SCR_MKT:  %f\n', SCR_MKT);
fprintf('SCR_LIFE: %f\n', SCR_LIFE);
fprintf('---------------------------------\n');

figure;
results = [BSCR(:)', SCR_MKT(:)', SCR_LIFE(:)'];
labels = {'BSCR', 'SCR\_MKT', 'SCR\_LIFE'};

bar(results)
set(gca, 'xticklabel', labels)
ylabel('Value')
title('Results')
text(1:numel(results), results, num2str(results', '%0.2f'), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom')

% Open the file for writing
fid = fopen('results.txt', 'a');

% Check if the file is opened successfully
if fid == -1
    error('Unable to open file for writing');
end

% Write the results to the file
fprintf(fid, 'Seed used: %d\n', Var_seed);
fprintf(fid, 'Number of simulations: %d\n', N);
fprintf(fid, 'BSCR: %f\n', BSCR);
fprintf(fid, 'SCR: %f\n', SCR);
fprintf(fid, 'SCR_MKT: %f\n', SCR_MKT);
fprintf(fid, 'SCR_LIFE: %f\n\n\n', SCR_LIFE);

% Close the file
fclose(fid);

disp('Results have been saved to "results.txt"');

% print results
fprintf('Results:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL);
fprintf('Death:  %f\n', Death_BEL);
fprintf('Expenses: %f\n', Expenses_BEL);
fprintf('Commission: %f\n', Commissions_BEL);
fprintf('---------------------------------\n');

fprintf('shock rates UP:\n');
fprintf('---------------------------------\n');

fprintf('Lapse:     %f\n', Lapse_BEL_rates_UP);
fprintf('Death:  %f\n', Death_BEL_rates_UP);
fprintf('Expenses: %f\n', Expenses_BEL_rates_UP);
fprintf('Commission: %f\n', Commissions_BEL_rates_UP);
fprintf('---------------------------------\n');


fprintf('shock rates DOWN:\n');
fprintf('---------------------------------\n');

fprintf('Lapse:     %f\n', Lapse_BEL_rates_DOWN);
fprintf('Death:  %f\n', Death_BEL_rates_DOWN);
fprintf('Expenses: %f\n', Expenses_BEL_rates_DOWN);
fprintf('Commission: %f\n', Commissions_BEL_rates_DOWN);
fprintf('---------------------------------\n');


fprintf('Equity risk:\n');
fprintf('---------------------------------\n');

fprintf('Lapse:     %f\n', Lapse_BEL_equity);
fprintf('Death:  %f\n', Death_BEL_equity);
fprintf('Expenses: %f\n', Expenses_BEL_equity);
fprintf('Commission: %f\n', Commissions_BEL_equity);
fprintf('---------------------------------\n');


fprintf('Property risk:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_property);
fprintf('Death:  %f\n', Death_BEL_property);
fprintf('Expenses: %f\n', Expenses_BEL_property);
fprintf('Commission: %f\n', Commissions_BEL_property);
fprintf('---------------------------------\n');


fprintf('Mortality risk:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_mortality);
fprintf('Death:  %f\n', Death_BEL_mortality);
fprintf('Expenses: %f\n', Expenses_BEL_mortality);
fprintf('Commission: %f\n', Commissions_BEL_mortality);
fprintf('---------------------------------\n');


fprintf('shock lapse UP:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_lapse_UP);
fprintf('Death:  %f\n', Death_BEL_lapse_UP);
fprintf('Expenses: %f\n', Expenses_BEL_lapse_UP);
fprintf('Commission: %f\n', Commissions_BEL_lapse_UP);
fprintf('---------------------------------\n');


fprintf('shock lapse DOWN:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_lapse_DOWN);
fprintf('Death:  %f\n', Death_BEL_lapse_DOWN);
fprintf('Expenses: %f\n', Expenses_BEL_lapse_DOWN);
fprintf('Commission: %f\n', Commissions_BEL_lapse_DOWN);
fprintf('---------------------------------\n');


fprintf('lapse mass risk:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_lapse_mass);
fprintf('Death:  %f\n', Death_BEL_lapse_mass);
fprintf('Expenses: %f\n', Expenses_BEL_lapse_mass);
fprintf('Commission: %f\n', Commissions_BEL_lapse_mass);
fprintf('---------------------------------\n');


fprintf('Expense risk:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_expense);
fprintf('Death:  %f\n', Death_BEL_expense);
fprintf('Expenses: %f\n', Expenses_BEL_expense);
fprintf('Commission: %f\n', Commissions_BEL_expense);
fprintf('---------------------------------\n');


fprintf('Catastrophe risk:\n');
fprintf('---------------------------------\n');
fprintf('Lapse:     %f\n', Lapse_BEL_cat);
fprintf('Death:  %f\n', Death_BEL_cat);
fprintf('Expenses: %f\n', Expenses_BEL_cat);
fprintf('Commission: %f\n', Commissions_BEL_cat);
fprintf('---------------------------------\n');