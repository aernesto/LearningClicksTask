% This script does the following:
%   1. generate the stimulus for 1 trial of the dynamic clicks task
%   2. evolve the full system with truncated gamma (CP count)
%   3. compute posterior over change point count at several time points
%   4. display a plot with two subplots described below
%       - raster plot of stimulus clicks, with change points marked
%       - posterior mean over CP count with 'error' bars shown at 1stdev
%
% REQUIRED SCRIPTS:
% returnPostA.m
% genStimBank2.m

clear
%% set parameters
% recall that h is normalized to 1 
hh=4;
snr=1;
% click rates in Hz
rateLow=0.01;
rateHigh=getlambdahigh(rateLow, snr, true);
% time step for forward Euler, in sec
dt=1/10000;
% max allowed change point count
gamma_max=100;
% hyperparameters for Gamma dist over hazard rate
alpha=1;
beta=1;
priorState=[.5,.5];

%trial duration (sec)
T=0.010; %10 msec

%% generate stimulus
% %data=genStimBank2(1,hh,rateHigh,rateLow,T);
% load('../data/ClickTrains_h4_rateHigh38_rateLow001_nTrials1_SHORT.mat')
% % number of distinct trial durations in the data array
% N=size(data,1);
% % get single trial from longest trial duration (3 sec)
% clicksCell=data{N,1};
% % total number of available trials for this trial duration
% nTrials = length(clicksCell);
% trial = 1; % select first trial for now
% [lTrain,rTrain, cptimes]=clicksCell{trial, 1:3};
lTrain=[];
rTrain=0.0049; % right before 5 msec
cptimes=0.004;


%% call ODE function
%posttimes=1:50;
posttimes=0.001:0.001:T;
%posttimes=[0.1:0.1:.9, posttimes];
posttimes(end)=T-2*dt;
% UP TO HEAR CODE IS FINE
tic
postH=returnPostH(lTrain, rTrain, rateLow, rateHigh, T, ...
    gamma_max, posttimes, priorState, alpha, beta, dt, cptimes);
toc



fig=figure(1); 
hax=axes; 
hold on 
SP=rTrain(1)*1000; %right click time in msec
%for snr=[.5,1,2,4]
%    rate_high=getlambdahigh(rate_low, snr, true);
%    P=jointPosteriorClicks(lTrain,rTrain);

    plot(1000*(posttimes), postH(1,:),'-o','LineWidth', ...
        3, 'MarkerSize',4)
%end
title('posterior prob of H+ as fcn of time')
ylabel('posterior prob(H+)')
xlabel('time within stimulus (msec)')
xlim([0,11])
ylim([0.45,1.05])
line([SP SP],get(hax,'YLim'),'Color',[1 0 0], 'LineWidth',2)
line(get(hax,'XLim'),[0.5,.5],'Color',[0 0 0], 'LineWidth',1)
line(get(hax,'XLim'),[1,1],'Color',[0 0 0], 'LineWidth',1)
legend('snr=1','click time', 'Location', 'east')
hax.FontSize=20;


% append prior values for time point t=0
% means=(0:gamma_max-1)*postGamma; % row vector of posterior means over gamma
% interm=((0:gamma_max-1).^2)*postGamma;
% stdevs=sqrt(interm-(means.^2));
% means=[0,means];
% stdevs=[0,stdevs];
% posttimes=[0,posttimes];
% %plot stimulus with change point times
% subplot(5,1,1)
% plotClickTrains(lTrain,rTrain,cptimes, rateLow, rateHigh)
% xlim([0,T])
% title('stimulus with CP times','FontSize',14)
% %plot posterior mean over gamma
% subplot(5,1,2)
% errorbar(posttimes,means,2*stdevs);
% xlim([0,T])
% title('posterior mean over gamma','FontSize',14)
% %plot posterior pmf at t=0 and first positive reporting time
% subplot(5,1,3)
% massOn0=.99;
% priorGamma=[massOn0,ones(1,gamma_max-1)*(1-massOn0)/gamma_max];
% bar(priorGamma)
% title('posterior at t=0')
% ylabel('probability')
% subplot(5,1,4)
% bar(postGamma(:,1))
% ylabel('probability')
% title(['posterior at t=',num2str(posttimes(2),2)])
% subplot(5,1,5)
% bar(postGamma(:,end))
% xlabel('CP count')
% ylabel('probability')
% title(['posterior at t=',num2str(posttimes(end),2)])