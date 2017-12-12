% This script generates a performance curve as a function of interrogation
% time, for the dynamic clicks task ODE model which assumes the hazard rate
% known.
% The data used for each trial is loaded from a .mat cell array.

%---------------------
% cell array structure
%---------------------
% cell name = data
% cell dimensions = N x 2;
%   N is the number of interrogation time points
%   right column contains the interrogation time
%   left column contains a cell array described below:
%       cell dimensions: M x 4
%           M = number of trials
%           column 1= left click train (col vector of click times)
%           column 2= right click train (col vector if click times)
%           column 3= change point times (col vector)
%           column 4= environment states 
%                   (col vector of 0's (for H-) and 1's (for H+))
%----------------------
tic
% load cell array
load('../data/ClickTrains_h1_rateHigh38_rateLow2_nTrials10000.mat')

% number of interrogation time points
N=length(data);

% click rates and hazard rate
rateHigh=38;
rateLow=2;
h=1;

% array that stores performance for each interrogation time
perf_array=zeros(N,1);

% loop over interrogation times
parfor tidx=1:N
    [clicksCell, interrogationTime]=data{tidx,:};

    % loop over trial and keep track of correctness of answer
    nTrials = length(clicksCell);
    correct = 0;
    for trial=1:nTrials
        trueEndState=clicksCell{trial, 4}(end);
        
        % get answer from ODE model
        [lTrain,rTrain]=clicksCell{trial, 1:2};
        answer=exactODE(lTrain, rTrain, rateLow, rateHigh, ...
            interrogationTime, h);
        % draw random choice if evidence is at 0 at end of trial
        if answer == 0
            rng('shuffle')
            if rand < 0.5
                answer = -1;
            else
                answer = 1;
            end
        end
        if and(answer == -1, trueEndState == 0)
            correct = correct + 1;
        elseif and(answer == 1, trueEndState == 1)
            correct = correct + 1;
        end
    end
    perf_array(tidx,1) = correct / nTrials; % nTrials
end
plot(.2:.2:.4,perf_array(1:2),'-o','linewidth',1)
%ylim([.5,1])
xlim([0,.4])
ylabel('performance', 'FontSize',20)
xlabel('interrogation time', 'FontSize',20)
%title('Performance of ODE h=1, rateHigh=38, rateLow=2','FontSize',20)
toc
        