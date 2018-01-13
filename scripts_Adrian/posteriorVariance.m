% This script does the following:
%   1. load the stimulus for 1 trial of the dynamic clicks task
%   2. evolve the full system with truncated gamma (CP count)
%   3. compute posterior variance over hazard rate at several time points

%% load stimulus
clear
load('../data/ClickTrains_h1_rateHigh38_rateLow2_nTrials1_LONG.mat')
% number of distinct trial durations in the data array
N=size(data,1);
% get single trial from longest trial duration (3 sec)
[clicksCell, T]=data{N,:};
%T
% total number of available trials for this trial duration
nTrials = length(clicksCell);
trial = 1; % select first trial for now
[lTrain,rTrain]=clicksCell{trial, 1:2};

%% set parameters
% recall that h is normalized to 1 and T was set above (3 sec)
% click rates in Hz
rateHigh=38;
rateLow=2;
% time step for forward Euler, in sec
dt=1/10000;
% max allowed change point count
gamma_max=100;
% hyperparameters for Gamma dist over hazard rate
alpha=1;
beta=1;
priorState=[.5,.5];
%% call ODE function
posttimes=1:50;
posttimes(end)=50-2*dt;
% UP TO HEAR CODE IS FINE
tic
[vars, means]=systemODE(lTrain, rTrain, rateLow, rateHigh, T, gamma_max,...
posttimes, priorState, alpha, beta, dt);
vars(vars==inf)=vars(1);
means(means==inf)=means(1);
toc
subplot(2,1,1)
plot(posttimes,means,'LineWidth',3)
ylabel('posterior mean','FontSize',14)
xlim([0,50])
ylim([0,2.2])
title('posterior mean over h','FontSize',14)
subplot(2,1,2)
plot(posttimes,vars,'LineWidth',3);
xlabel('time','FontSize',14)
ylabel('posterior var','FontSize',14)
xlim([0,50])
ylim([0,max(abs(vars))]);
title('posterior var over h','FontSize',14)