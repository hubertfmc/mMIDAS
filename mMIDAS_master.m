
clear
clc
%%
[DataY,DataYdate] = xlsread('mydata.xlsx','sheet1');
DataYdate = DataYdate(2:end,1);
[DataX,DataXdate] = xlsread('mydata.xlsx','sheet2');
DataXdate = DataXdate(2:end,1);

DataXgrowth = log(DataX(2:end)./DataX(1:end-1))*100;
DataYgrowth = log(DataY(2:end)./DataY(1:end-1))*100;
DataX=DataXgrowth;
DataY=DataYgrowth;
DataYdate = DataYdate(2:end);
DataXdate = DataXdate(2:end);

% Specify lag structure and sample size 
Xlag = 9;
Ylag = 1;
% Horizon = 1;
Horizon = 3;
EstStart = '1985-01-01';
EstEnd = '2009-01-01';
Method = 'fixedWindow';
%% (1) Univariate setting
[OutputForecast3,OutputEstimate3] = m1MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'Display','estimate');
%% (2) Multivariate setting with hypothetical data
% Adding random data for testing
nx = 1; % Specify the number of new variables in matrix X;
DataX = [DataX, randn(size(DataX,1),nx)];
Ylag = 3; % Including more lag Y
[OutputForecast3,OutputEstimate3] = m1MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,...
    'Xlag',Xlag,'Ylag',Ylag,'Horizon',Horizon,'EstStart',EstStart,'EstEnd',EstEnd,'Polynomial','expAlmon','Method',Method,'Display','estimate');
