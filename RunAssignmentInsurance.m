% initial settings
clc
clear all
close all

% fix random seed
rng(42)

% start run time
tic

% load data from xls file

% filename in .xls
filename = 'LifeTable2022_2023.xls';

% read excel data from filename
% INPUT filename, formatData, interval of rows and columns to read e.g. 'E8:F36'
% 
survivors = xlsread(filename,1,'D8:D36');
deaths = xlsread(filename,1,'F8:F31');
P_death = xlsread(filename,1,'H8:H31');
Years_lived = xlsread(filename,1,'J8:J31');
Proj_prob = xlsread(filename,1,'L8:L31');
Life_exp = xlsread(filename,1,'N8:N31');

filename = 'EIOPA_RFR_20240331_Term_Structures.xlsx';

% disp(survivors)
% disp(deaths)
% disp(P_death)
% disp(Years_lived)
% disp(Proj_prob)
% disp(Life_exp)


Sheet = 'RFR_spot_no_VA';

rates = xlsread(filename,1,'S11:S160');

disp(rates)

% plot yield rate curve
figure
years = 1:150;
plot(years,rates)

% plot mortality table
figure
plot([1:24],Life_exp)

%tolo gay
