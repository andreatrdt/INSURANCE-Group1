%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment Insurance
% Authors:
%   - Giacomo Manfredi
%   - Lorenzo Tolomelli
%   - Andera Tarditi
%   - Giovanni Riondato
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial settings
clc;
clear all;
close all;

% select format
format bank

% fix random seed
Var_seed = 0; % the answer to the ultimate question of life, the universe, and everything
rng(Var_seed)

% start run time
tic

% load data from xlsx file

% filename in .xls
filename = 'LifeTable.xlsx';

% read excel data from LifeTable.xlsx
P_death = xlsread(filename,1,'D64:D113');
P_death = P_death./1000;


load('prob_death.mat');
P_death = prob_death;


% read excel data from EIOPA
filename = 'EIOPA_RFR_20240331_Term_Structures.xlsx';

rates_UP = xlsread(filename,5,'S11:S160');
rates_DOWN = xlsread(filename,6,'S11:S160');
rates = xlsread(filename,1,'S11:S160');

% convert rates to log returns
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

% parameters

F0 = 1e5; % value of fund at t = 0
S0 = 0.8*F0; % value of equity at t = 0
sigma_equity = 0.2; % volatility of equity
T = 50; % number of years
N = 1e6; % number of simulations
regular_deduction = 0.022; % regular deduction
sigma_pf = 0.1; % volatility of property features
yearly_expense_t0 = 50; % yearly expenses at t = 0
dt = [1:1:T]';  % time grid
inflation = 0.02; % inflation rate
expenses = yearly_expense_t0.*(1+inflation).^[0; dt(1:end-1)]; % expenses
benefit_commission = 20; % benefit commission

disp('Parameters loaded...')

% Liabilities setup
COMM = 0.014; % commission
lt = 0.15 * ones(length(dt),1); % lapse rate

% consider only the first T rates
rates = rates(1:T);
rates_UP = rates_UP(1:T);
rates_DOWN = rates_DOWN(1:T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intrest rate risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Interest rate risk analysis...')

%% Basic scenario

% simulate equity prices
S_simulated = simulate_GBM(rates, S0, sigma_equity, T, N, regular_deduction);
S = mean(S_simulated,1);


% simulate property features
PF_0 = 0.2*F0;
PF_simulated = simulate_GBM(rates, PF_0, sigma_pf, T, N, regular_deduction);
PF = mean(PF_simulated,1);

% calculate fund value
F = S + PF;

%calculate discouts
discounts = exp(-rates.*dt);

liabilities = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses, dt, F, benefit_commission, T);

BOF = F0 - liabilities;

disp(BOF)

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
S_rates_UP = mean(S_rates_UP,1);
% simulate property features
PF_rates_UP = simulate_GBM(rates_UP, PF_0, sigma_pf, T, N, regular_deduction);
PF_rates_UP = mean(PF_rates_UP, 1);

% calculate fund value
F_rates_UP = [F0 , S_rates_UP + PF_rates_UP];

%calculate discouts
discounts_UP = exp(-rates_UP.*dt);

% Liabilities
liabilities_rates_UP = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt, F_rates_UP, benefit_commission, T);

BOF_rates_UP = F0 - liabilities_rates_UP;
delta_BOF_rates_UP = max(BOF-BOF_rates_UP,0);

%% Stress scenario DOWN

% simulate equity prices
S_rates_DOWN = simulate_GBM(rates_DOWN, S0, sigma_equity, T, N, regular_deduction);
S_rates_DOWN = mean(S_rates_DOWN,1);
% simulate property features
PF_rates_DOWN = simulate_GBM(rates_DOWN, PF_0, sigma_pf, T, N, regular_deduction);
PF_rates_DOWN = mean(PF_rates_DOWN,1);

%calculate discouts
discounts_DOWN = exp(-rates_DOWN.*dt);

% calculate fund value
F_rates_DOWN = S_rates_DOWN + PF_rates_DOWN;

% Liabilities
liabilities_rates_DOWN = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F_rates_DOWN, benefit_commission,T);

BOF_rates_DOWN = F0 - liabilities_rates_DOWN;
delta_BOF_rates_DOWN = max(BOF-BOF_rates_DOWN,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Property risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Property risk analysis...')

% Initial value of the shocked value of the property
P0_shocked=(1-0.25)*PF_0;
% initial value for the fund
F0_property = P0_shocked + S0;  

% Property simulation
P_shocked = simulate_GBM(rates(1:T), P0_shocked, sigma_pf, T, N, regular_deduction);
P_shocked = mean(P_shocked,1);

% Value of the fund at each time step
F = P_shocked + S;   

% Computation of Liabilities 
liabilities_shocked_pr = Liabilities(F0_property, P_death, lt, regular_deduction, COMM, discounts, expenses, dt, F, benefit_commission, T);

% Delta BOF
BOF_pr = F0_property - liabilities_shocked_pr;
delta_BOF_pr = max(BOF - BOF_pr,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equity risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Equity risk analysis...')

% initial value for shocked equity
S0_shocked = (1-0.39)*S0;
% initail value for the fund
F0_equity = S0_shocked + PF_0;

% Equity simulation
S_shocked = simulate_GBM(rates(1:T), S0_shocked, sigma_equity, T, N, regular_deduction);
S_shocked = mean(S_shocked,1);

F = S_shocked + PF;       % new value of the fund at each time step

% Computation of Liabilities 
liabilities_shocked_eq = Liabilities(F0_equity, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% delta BOF
BOF_eq = F0_equity - liabilities_shocked_eq;
delta_BOF_eq = max(BOF - BOF_eq,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORTALITY RISK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Mortality risk analysis...')

P_death_shocked = min(1,P_death*1.15);

% simulate equity prices
F = S + PF;

% Computation of Liabilities
liabilities_shocked_mortality = Liabilities(F0, P_death_shocked, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
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

% Computation of Liabilities
liabilities_shocked_lapse_UP = Liabilities(F0, P_death, lt_shocked_UP, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
BOF_lapse_UP = F0 - liabilities_shocked_lapse_UP;
delta_BOF_lapse_UP = max(BOF - BOF_lapse_UP,0);

%% lapse DOWN risk

lt_shocked_DOWN = max(0.5*lt,lt-0.2) * ones(length(dt),1);

% Computation of Liabilities
liabilities_shocked_lapse_DOWN = Liabilities(F0, P_death, lt_shocked_DOWN, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
BOF_lapse_DOWN = F0 - liabilities_shocked_lapse_DOWN;
delta_BOF_lapse_DOWN = max(BOF - BOF_lapse_DOWN,0);

%% Lapse mass risk

lt_shocked_mass = [lt+0.4;lt*ones(T-1,1)];

% Computation of Liabilities
liabilities_shocked_lapse_mass = Liabilities(F0, P_death, lt_shocked_mass, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
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

% Computation of Liabilities
liabilities_shocked_expense = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses_shocked,dt,F,benefit_commission,T);

% Delta BOF
BOF_expense = F0 - liabilities_shocked_expense;
delta_BOF_expense = max(BOF - BOF_expense,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Catastrophe risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Catastrophe risk analysis...')

% Catastrophe scenario
P_death_cat = [P_death(1)+0.0015; P_death(2:end)];

% Computation of Liabilities 
liabilities_cat = Liabilities(F0, P_death_cat, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

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
SCR_IR = max(delta_BOF_rates_UP , delta_BOF_rates_DOWN);        % Interest rate risk
SCR_EQ = delta_BOF_eq;                              % Equity risk
SCR_PR = delta_BOF_pr;                              % Property risk 
MKT_vec = [SCR_IR ; SCR_EQ ; SCR_PR]; 

% check if we are exposed to IR up or down to choose the correlation matrix
if delta_BOF_rates_UP > delta_BOF_rates_DOWN
    MKT_corr = [1 0 0; 0 1 0.75; 0 0.75 1];
else
    MKT_corr = [1 0.5 0.5; 0.5 1 0.75; 0.5 0.75 1];
end
SCR_MKT = sqrt(MKT_vec' * MKT_corr * MKT_vec);

% SCR Life risk
SCR_MORT = delta_BOF_mortality;                                                     % Mortality risk
SCR_LAPSE = max([delta_BOF_lapse_UP,delta_BOF_lapse_DOWN,delta_BOF_lapse_mass]);    % Lapse risk
SCR_EXP = delta_BOF_expense;                                                        % Expense risk
SCR_CAT = delta_BOF_catastrophe;                                                    % Catastrophe risk
LIFE_vec = [SCR_MORT ; SCR_LAPSE ; SCR_EXP ; SCR_CAT];

% correlation matrix
LIFE_corr = [1 0 0.25 0.25; 0 1 0.5 0.25; 0.25 0.5 1 0.25 ; 0.25 0.25 0.25 1 ];

SCR_LIFE = sqrt(LIFE_vec' * LIFE_corr * LIFE_vec);

% BSCR
SCR = [SCR_MKT ; SCR_LIFE];
SCR_corr = [1 0.25; 0.25 1];        % correlation matrix
BSCR = sqrt(SCR' * SCR_corr * SCR);

% print results

fprintf('Results:\n');
fprintf('---------------------------------\n');
fprintf('BSCR:     %f\n', BSCR);
fprintf('SCR_MKT:  %f\n', SCR_MKT);
fprintf('SCR_LIFE: %f\n', SCR_LIFE);
fprintf('---------------------------------\n');

figure;
results = [BSCR(:)', SCR(:)', SCR_MKT(:)', SCR_LIFE(:)'];
labels = {'BSCR', 'SCR', 'SCR\_MKT', 'SCR\_LIFE'};

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

% end run time
toc