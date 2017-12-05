function [ct,E]=genClickEnvt(stimulusLength, h)
%arguments: 
    %stimulusLength is the stimulus length in seconds
    %h is the hazard rate in Hz
%returns two column vectors: 
    %ct contains the change point times
    %E of length length(ct)+1, contains the states H+ and H-, as 1 and 0
    %respectively.
%IMPORTANT: this function draws the initial environment state uniformly.
%That is, H+ and H- have equal probability of happening at time 0.

ct=genPoissonTrain(h, stimulusLength); % generate change point times

E=zeros(length(ct)+1,1);
idx=1; % initialize idx for coming while loop

% draw initial state uniformly
initialState=round(rand); % 0 codes for H- ,1 for H+
E(idx)=initialState;

while idx < length(E)
    E(idx+1)=1-E(idx);  % switch state after each change point
    idx=idx+1;
end
end