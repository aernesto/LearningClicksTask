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

fig=figure(1); 
SP=rTrain(1)*1000; %right click time in msec
ax1=subplot(1,2,1);
hold on
title('post H+')
ylabel('prob(H+)')
xlabel('msec')
xlim([0,11])
ylim([0,1.05])
line([SP SP],get(ax1,'YLim'),'Color',[1 0 0], 'LineWidth',2)
line(get(ax1,'XLim'),[0.5,.5],'Color',[0 0 0], 'LineWidth',1)
line(get(ax1,'XLim'),[1,1],'Color',[0 0 0], 'LineWidth',1)
ax1.FontSize=20;

ax2=subplot(1,2,2);
hold on 
title('post H-')
ylabel('prob(H-)')
xlabel('msec')
xlim([0,11])
ylim([0,1.05])
line([SP SP],get(ax2,'YLim'),'Color',[1 0 0], 'LineWidth',2)
line(get(ax2,'XLim'),[0.5,.5],'Color',[0 0 0], 'LineWidth',1)
line(get(ax2,'XLim'),[1,1],'Color',[0 0 0], 'LineWidth',1)
ax2.FontSize=20;

for snr=[.5,1,2,4]
    rateHigh=getlambdahigh(rateLow, snr, true);
    postH=returnPostH(lTrain, rTrain, rateLow, rateHigh, T, ...
    gamma_max, posttimes, priorState, alpha, beta, dt, cptimes);
    plot(ax1,1000*posttimes, postH(1,:), '-*','LineWidth',3)
    plot(ax2,1000*posttimes, postH(2,:), '-*','LineWidth',3)
%     ax1=subplot(2,2,1);
%     ax1.FontSize=20;
%     theoPost=1./(posttimes+1);
%     CCC1=jointPost(1,1)/theoPost(1);
%     CCC2=jointPost(1,5)/theoPost(5);
%     CCC3=jointPost(gamma_max+1,5)/theoPost(5);
%     theoPost1=[CCC1*theoPost(1:4),CCC2*theoPost(5:end)];
%     theoPost2=[CCC1*theoPost(1:4),CCC3*theoPost(5:end)];
%     hold on
%     plot(1000*(posttimes), jointPost(1,:),'bo',...
%         1000*posttimes, theoPost1, 'r*','MarkerSize',10)
%     line([SP SP],get(ax1,'YLim'),'Color',[1 0 0], 'LineWidth',2)
%     line(get(ax1,'XLim'),[0.5,.5],'Color',[0 0 0], 'LineWidth',1)
%     line(get(ax1,'XLim'),[1,1],'Color',[0 0 0], 'LineWidth',1)
%     legend('sim','theo')
%     title('H+, a=0')
%     ax2=subplot(2,2,2);
%     ax2.FontSize=20;
%     hold on
%     plot(1000*(posttimes), jointPost(gamma_max+1,:),'bo',...
%         1000*posttimes, theoPost2, 'r*','MarkerSize',10)
%     %legend('sim','theo','Location','northeast')
%     line([SP SP],[0,1],'Color',[1 0 0], 'LineWidth',2)
%     line(get(ax2,'XLim'),[0.5,.5],'Color',[0 0 0], 'LineWidth',1)
%     line(get(ax2,'XLim'),[1,1],'Color',[0 0 0], 'LineWidth',1)
%     title('H-, a=0')
%     ax3=subplot(2,2,3);
%     ax3.FontSize=20;
%     plot(1000*(posttimes), jointPost(2,:),'-o','LineWidth', ...
%         3, 'MarkerSize',8)
%     title('H+, a=1')
%     ax4=subplot(2,2,4);
%     ax4.FontSize=20;
%     plot(1000*(posttimes), jointPost(gamma_max+2,:),'-o','LineWidth', ...
%         3, 'MarkerSize',8)
%     title('H-, a=1')

end


legend(ax1,'click time','','','snr=0.5','snr=1','snr=2','snr=4', 'Location', 'east')
legend(ax2,'click time','','','snr=0.5','snr=1','snr=2','snr=4', 'Location', 'east')



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