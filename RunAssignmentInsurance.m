% initial settings
clc;
clear all;
close all;

% fix random seed
rng(42)

% start run time
tic

% load data from xls file

% filename in .xls
filename = 'LifeTable.xlsx';

% read excel data from filename
% INPUT filename, formatData, interval of rows and columns to read e.g. 'E8:F36'
% 
P_death = xlsread(filename,1,'D64:D113');
P_death = P_death./1000;

filename = 'EIOPA_RFR_20240331_Term_Structures.xlsx';

Sheet = 'RFR_spot_no_VA';

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

F0=1e5; % value of the funf at t = 0

% load data
S0=0.8*F0;
sigma_equity=0.2;
T = 50;
N = 10;
regular_deduction = 0.022;
sigma_pf = 0.1;
yearly_expense_t0 = 50;
dt = [1:1:T]'; 
inflation = 0.02;
expenses = yearly_expense_t0.*(1+inflation).^[0; dt(1:end-1)];
rates = rates(1:T);
rates_UP = rates_UP(1:T);
rates_DOWN = rates_DOWN(1:T);
benefit_commission = 20;

% Liabilities
COMM = 0.014;
lt = 0.15 * ones(length(dt),1);

%% Basic scenario

S = simulate_GBM(rates(1:T), S0, sigma_equity, T, N, regular_deduction);
S=mean(S,1);
% simulate property features
PF_0 = 0.2*F0;
PF = simulate_GBM(rates(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF = mean(PF,1);

F= S + PF;

%calculate discouts
discounts = exp(-rates.*dt);

liabilities = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F,benefit_commission,T);
disp(liabilities)

%% Stress scenario UP

%simulate equity prices
S_UP = simulate_GBM(rates_UP(1:T), S0, sigma_equity, T, N, regular_deduction);
S_UP = mean(S_UP,1);
PF_UP = simulate_GBM(rates_UP(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF_UP = mean(PF_UP,1);

F_UP = S_UP + PF_UP;

%calculate discouts
discounts_UP = exp(-rates_UP.*dt);


% Liabilities
liabilities_UP = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F_UP,benefit_commission,T);
disp(liabilities_UP)


%% Stress scenario DOWN

S_DOWN = simulate_GBM(rates_DOWN(1:T), S0, sigma_equity, T, N, regular_deduction);
S_DOWN = mean(S_DOWN,1);
PF_DOWN = simulate_GBM(rates_DOWN(1:T), PF_0, sigma_pf, T, N, regular_deduction);
PF_DOWN = mean(PF_DOWN,1);

%calculate discouts
discounts_DOWN = exp(-rates_DOWN.*dt);

F_DOWN = S_DOWN + PF_DOWN;

% Liabilities
liabilities_DOWN = Liabilities(F0, P_death, lt, regular_deduction, COMM, discounts, expenses,dt,F_DOWN, benefit_commission,T);
disp(liabilities_DOWN)

toc

