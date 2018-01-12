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
[lTrain,rTrain]=clicksCell{trial, 1:2};

%% set parameters
% recall that h is normalized to 1 and T was set above (3 sec)
% click rates in Hz
rateHigh=38;
rateLow=2;
% time step for forward Euler, in sec
dt=1/10000;
% max allowed change point count
gamma_max=50;
% hyperparameters for Gamma dist over hazard rate
alpha=1;
beta=1;
priorState=[.5,.5];
%% call ODE function
posttimes=.5:.5:3;
posttimes(end)=3-2*dt;
tic
vars=systemODE(lTrain, rTrain, rateLow, rateHigh, T, gamma_max,...
posttimes, priorState, alpha, beta, dt);
toc
plot(posttimes,vars);
xlabel('time')
ylabel('posterior var')
ylim([0,max(vars)]);
title('posterior variance over h single trial')