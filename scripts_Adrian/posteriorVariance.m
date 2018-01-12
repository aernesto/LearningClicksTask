% This script does the following:
%   1. load the stimulus for 1 trial of the dynamic clicks task
%   2. evolve the full system with truncated gamma (CP count)
%   3. compute posterior variance over hazard rate at several time points

%% load stimulus
clear
load('../data/ClickTrains_h1_rateHigh38_rateLow2_nTrials10000.mat')
% number of distinct trial durations in the data array
N=length(data);
% get single trial from longest trial duration (3 sec)
[clicksCell, T]=data{N,:};
% total number of available trials for this trial duration
nTrials = length(clicksCell);
trial = 1; % select first trial for now
[lTrain1,rTrain1]=clicksCell{1, 1:2};

[lTrain2,rTrain2]=clicksCell{2, 1:2};
[lTrain3,rTrain3]=clicksCell{3, 1:2};
[lTrain4,rTrain4]=clicksCell{4, 1:2};
[lTrain5,rTrain5]=clicksCell{5, 1:2};
[lTrain6,rTrain6]=clicksCell{6, 1:2};
lTrain=[lTrain1;3+lTrain2;6+lTrain3;9+lTrain4;12+lTrain5;15+lTrain6];
rTrain=[rTrain1;3+rTrain2;6+rTrain3;9+rTrain4;12+rTrain5;15+rTrain6];

%% set parameters
% recall that h is normalized to 1 and T was set above (3 sec)
% click rates in Hz
rateHigh=38;
rateLow=2;
% time step for forward Euler, in sec
dt=1/10000;
% max allowed change point count
gamma_max=150;
% hyperparameters for Gamma dist over hazard rate
alpha=1;
beta=1;
priorState=[.5,.5];
%% call ODE function
posttimes=.5:.5:18;
posttimes(end)=18-2*dt;
T=18;
tic
vars=systemODE(lTrain, rTrain, rateLow, rateHigh, T, gamma_max,...
posttimes, priorState, alpha, beta, dt);
toc
plot(posttimes,vars);
xlabel('time')
ylabel('posterior var')
ylim([0,max(vars)]);
title('posterior variance over h single trial')