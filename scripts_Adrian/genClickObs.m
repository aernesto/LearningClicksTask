function [lTrain,rTrain]=genClickObs(ct,E, rateLow, rateHigh, stimulusLength)
% generates the two trains of clicks that the rat hears
%arguments:
%   ct: is a column vector with change point times
%   E: is a column vector of length length(ct)+1 with the state H+/- of the
%   environment encoded as binary value (1 for H+)
%   rateLow: low click rate in Hz
%   rateHigh: high click rate in Hz
%   stimulusLength: stimulus length in sec
% returns: two column vectors for left and right trains respectively

% We use the convention, H+ corresponds to rateHigh to the right ear
% We use the memoryless property of the Poisson trains to simply 'stack'
% them, after each change point.

if ~((length(ct) + 1) == length(E))
    error('change point times vector should have one more entry than environment vector')
end

nTrains=length(E); % number of trains to stack, for each ear

trains=cell(nTrains,2); % cell array storing event trains for each ear
                        % column 1/2 for train left/right
                        
%construct trains between each change point
for tt=1:nTrains
    % extract time length of current train
    
    % 'timeLength' is the length of the time interval on which the click rate
    % is constant
    % 'offset' is the time, in seconds, at which the interval starts
    
    % case: no change point in whole trial
    if nTrains == 1
        timeLength = stimulusLength;
        offset = 0;
    % case: first interval in a trial with more than 0 change points
    elseif tt == 1
        timeLength = ct(tt);
        offset = 0;
    % case: last interval in a trial with more than 1 interval
    elseif tt == nTrains
        offset = ct(end);
        timeLength = stimulusLength - offset;    
    % case: intermediate interval, in a trial with more than 2 intervals
    else
        offset = ct(tt-1);
        timeLength = ct(tt) - offset;    
    end
        
    % construct trains for both ears, depending on envt state
    
    % case: envt is in state H+ --> high rate to right ear
    if E(tt) 
        trains{tt,2}=genPoissonTrain(rateHigh,timeLength)+ offset; % right ear
        trains{tt,1}=genPoissonTrain(rateLow,timeLength) + offset;  % left ear
    % case: envt in state H- ---> high rate to left ear
    else 
        trains{tt,2}=genPoissonTrain(rateLow,timeLength) + offset;  % right ear
        trains{tt,1}=genPoissonTrain(rateHigh,timeLength) + offset; % left ear
    end
end

% concatenate trains. Each train is a column vector
lTrain = cell2mat(trains(:,1));
rTrain = cell2mat(trains(:,2));
end