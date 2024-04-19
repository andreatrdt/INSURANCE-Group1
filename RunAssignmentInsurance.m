%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 2: Insurance
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

% select format
format bank

% fix random seed
rng(42)

% start run time
tic

% load data from xls file

% filename in .xls
filename = 'LifeTable.xlsx';

% read excel data from LifeTable.xlsx
P_death = xlsread(filename,1,'D64:D113');
P_death = P_death./1000;

% read excel data from EIOPA
filename = 'EIOPA_RFR_20240331_Term_Structures.xlsx';

rates_UP = xlsread(filename,5,'S11:S160');
rates_DOWN = xlsread(filename,6,'S11:S160');
rates = xlsread(filename,1,'S11:S160');


% plot yield rate curve
figure
hold on 
years = 1:150;
plot(years,rates)
plot(years,rates_UP)
plot(years,rates_DOWN)
rate_50 = interp1(years,rates,50);
plot(50,rate_50,'ro')
legend( 'rates','rates_{DOWN}','rates_{UP}')

% parameters

F0 = 1e5; % value of fund at t = 0
S0 = 0.8*F0; % value of equity at t = 0
sigma_equity = 0.2; % volatility of equity
T = 50; % number of years
N = 1e4; % number of simulations
regular_deduction = 0.022; % regular deduction
sigma_pf = 0.1; % volatility of property features
yearly_expense_t0 = 50; % yearly expenses at t = 0
dt = [1:1:T]';  % time grid
inflation = 0.02; % inflation rate
expenses = yearly_expense_t0.*(1+inflation).^[0; dt(1:end-1)]; % expenses
benefit_commission = 20; % benefit commission

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

%% Basic scenario

% simulate equity prices
S = simulate_GBM(rates(1:T), S0, sigma_equity, T, N, regular_deduction);
S = mean(S,1);
% simulate property features
PF_0 = 0.2*F0;
PF = simulate_GBM(rates(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF = mean(PF,1);

% calculate fund value
F = S + PF;

%calculate discouts
discounts = exp(-rates.*dt);

liabilities = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);
disp(liabilities)

BOF = F0 - liabilities;

%% Stress scenario UP

%simulate equity prices
S_UP = simulate_GBM(rates_UP(1:T), S0, sigma_equity, T, N, regular_deduction);
S_UP = mean(S_UP,1);
% simulate property features
PF_UP = simulate_GBM(rates_UP(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF_UP = mean(PF_UP,1);

% calculate fund value
F_UP = S_UP + PF_UP;

%calculate discouts
discounts_UP = exp(-rates_UP.*dt);

% Liabilities
liabilities_UP = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F_UP,benefit_commission,T);
disp(liabilities_UP)

BOF_UP = F0 - liabilities_UP;
delta_BOF_UP = max(BOF-BOF_UP,0);

%% Stress scenario DOWN

% simulate equity prices
S_DOWN = simulate_GBM(rates_DOWN(1:T), S0, sigma_equity, T, N, regular_deduction);
S_DOWN = mean(S_DOWN,1);
% simulate property features
PF_DOWN = simulate_GBM(rates_DOWN(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF_DOWN = mean(PF_DOWN,1);

%calculate discouts
discounts_DOWN = exp(-rates_DOWN.*dt);

% calculate fund value
F_DOWN = S_DOWN + PF_DOWN;

% Liabilities
liabilities_DOWN = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F_DOWN, benefit_commission,T);
disp(liabilities_DOWN)

BOF_DOWN = F0 - liabilities_DOWN;
delta_BOF_DOWN = max(BOF-BOF_DOWN,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Property risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial value of the shocked value of the property
P0_shocked=(1-0.25)*PF_0;
% initial value for the fund
F0_property = P0_shocked + S0;    

% Property simulation
P_shocked = simulate_GBM(rates(1:T), P0_shocked, sigma_pf, T, N, regular_deduction);

% Value of the fund at each time step
F = P_shocked + S;       

% Computation of Liabilities 
liabilities_shocked_pr = Liabilities(F0_property, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
BOF_pr = F0_property - liabilities_shocked_pr;
delta_BOF_pr = max(BOF - BOF_pr,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equity risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial value for shocked equity
S0_shocked = (1-0.39)*S0;
% initail value for the fund
F0_equity = S0_shocked + PF_0;

% Equity simulation
S_shocked = simulate_GBM(rates(1:T), S0_shocked, sigma_equity, T, N, regular_deduction);

F = S_shocked + PF;       % new value of the fund at each time step

% Computation of Liabilities 
liabilities_shocked_eq = Liabilities(F0_equity, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% delta BOF
BOF_eq = F0_equity - liabilities_shocked_eq;
delta_BOF_eq = max(BOF - BOF_eq,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORTALITY RISK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_death_shocked = min(1,P_death*1.15);

% Computation of Liabilities
liabilities_shocked_mortality = Liabilities(F0, P_death_shocked, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);

% Delta BOF
BOF_mortality = F0 - liabilities_shocked_mortality;
delta_BOF_mortality = max(BOF - BOF_mortality,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lapse risk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% SCR Market risk
SCR_IR = max(delta_BOF_UP , delta_BOF_DOWN);        % Interest rate risk
SCR_EQ = delta_BOF_eq;                              % Equity risk
SCR_PR = delta_BOF_pr;                              % Property risk 
MKT_vec = [SCR_IR ; SCR_EQ ; SCR_PR]; 

% check if we are exposed to IR up or down to choose the correlation matrix
if delta_BOF_UP > delta_BOF_DOWN
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
fprintf('BSCR: %f\n', BSCR)
fprintf('SCR: %f\n', SCR)
fprintf('SCR_MKT: %f\n', SCR_MKT)
fprintf('SCR_LIFE: %f\n', SCR_LIFE)

% end run time
toc

