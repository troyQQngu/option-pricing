% MATH 512 proj 5
clear all
close all

%% Preparation
TSLA = readmatrix("TSLA.csv");
S = TSLA(:,6); % adjusted closing price

% basic parameter
T = 0.1587;
S0 = 1004.3;
K = S0 + 200;
r = 0.02;

% compute historical daily return
mu = zeros(length(S)-1,1);
for i = 1:length(mu)
    mu(i) = log(S(i+1)/S(i));
end

%compute historical daily volatility

sig_his = std(mu);
sig_an = sig_his*sqrt(252);
%% Question f
% monte carlo simulation 

L = 100; % number of partitions
n = 20000; % number of sample

m = 100;% to test variance
t=tic;
C0 = zeros(m,1);
for i = 1:m
    S = mcsim(S0,T,L,r,sig_an,n);
    payoff = max([S-K;zeros(1,n)],[],1);
    C0(i) = exp(-r*T)*mean(payoff);
end
rt_mc = toc(t);

std_mc = std(C0);

%% Question g
% Antithetic variates


t=tic;
C0_anti = zeros(m,1);
for i = 1:m
    S_anti = mcAntithetic(S0,r,sig_an,T,L,n);
    payoff_anti = max([S_anti-K;zeros(1,2*n)],[],1);
    C0_anti(i) = exp(-r*T)*mean(payoff_anti);
end
rt_av = toc(t);

std_av = std(C0_anti);
%% Question h
% Control Variates. We will use the built-in function blsprice as the
% 'accurate' solution


t=tic;
C0_cv = zeros(m,1);
for i = 1:m
    C0_cv(i) = mcControlVar(S0,K,r,T,sig_an,L,n);
end
rt_cv = toc(t);

std_cv = std(C0_cv);
%% Question i
% implied volatility

blsimpv(S0,K,r,T,26.56);

