function t=genPoissonTrain(rate,timeLength)
%generates poisson train of 'timeLength' seconds with rate 'rate' in Hz, 
%using the Gillespie algorithm.
%returns:
    %t is a column vector containing event times in second
    %%return empty vector if no event occurred in time interval

%pre-allocate 100 times the mean array size 
%for speed, will be shrinked after computation
t=zeros(ceil(100*rate*timeLength),1); % vector containing the event times

totalTime=0;
eventIndex=0;
while totalTime < timeLength
    sojournTime=exprnd(1/rate);
    totalTime=totalTime+sojournTime;
    eventIndex=eventIndex+1;
    t(eventIndex) = totalTime;
end
%-----
%NOTE: by now, t must have at least its first entry non zero
%-----

%trim unused nodes, and maybe last event if occurred beyond timeLength
[lastEvent,idxLastEvent]=max(t);
%-----
%NOTE: idxLastEvent is at least 1
%-----
if lastEvent > timeLength
    idxLastEvent = idxLastEvent - 1;
end

if idxLastEvent == 0
    t=zeros(0,1); %return empty vector if no event occurred in time interval
else
    t=t(1:idxLastEvent);
end
%t
end
