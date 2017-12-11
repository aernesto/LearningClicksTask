%This script uses parallel processing to generate and store a bank
%of independent trials from the dynamic clicks task
%returns: a file named ClickTrains.mat containing a single cell array
%           named data. This cell array has dimensions Nx2, where N is the
%           total number of trial durations for which trials should be
%           generated. The first column contains a cell array of trials. 
%           The second column contains the trial duration in seconds.
%           The aforementioned 'array of trials' has dimentions Mx4, where
%           M is the number of trials and 1st and 2st columns contain
%           the left and right click trains, respectively. Column 3
%           contains the change point times and column 4 contains the
%           environment state.

%parameters
nTrials = 5;    % number of trials for each trial duration
h = 1;          % hazard rate in Hz for environmental changes
rateHigh = 38;  % highest click rate in Hz
rateLow = 2;    % lowest click rate in Hz
interrogationTimes = 0.2:0.2:3; % vector of interrogation times
totIdx = length(interrogationTimes);
data = cell(totIdx,2);

for idx = 1:totIdx
    data{idx, 2} = interrogationTimes(idx);
end
parpool(16)

parfor idx = 1:totIdx
    T = interrogationTimes(idx);
    rng('shuffle')
    trains=cell(nTrials,4);
    for jj = 1:nTrials
        [ct,E]=genClickEnvt(T, h);
        trains{jj,3} = ct;
        trains{jj,4} = E;
        [lTrain, rTrain]=genClickObs(ct, E, rateLow, rateHigh, T);
        trains{jj,1} = lTrain;
        trains{jj,2} = rTrain;
    end
    data{idx,1} = trains;
end

save('/scratch/adrian/ClickTrains_h1_rateHigh38_rateLow2_nTrials5.mat',...
'data','-v7.3')
