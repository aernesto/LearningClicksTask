function ss=exactODE(lTrain, rTrain, rateLow, rateHigh, stimulusLength, h)
% This function evolves the jump ODE for the known hazard rate, over a 
% given stimulus train, and returns the sign of the evidence at the end
% of the trial

% maximum value of the indices for each train vector
lmax=length(lTrain);
rmax=length(rTrain);

% indices of the next entry to check in each train
ltidx=0;
rtidx=0;

% set times of next left and right clicks (infinity if non-existent)
if lmax > 0
    ltidx=1;
    lnxt=lTrain(1);
else
    lnxt=inf;
end
if rmax > 0
    rtidx=1;
    rnxt=rTrain(1);
else
    rnxt=inf;
end

% initialize evidence amount and time
evidence=0;
time=0;

% set delta jump sizes
kappa = log(rateHigh/rateLow);

% to evolve the ODE, the algorithm from time t up to t_new is as follows:
% set t_new as the minimum of the three following time:
    % time of next right click (set to infinity if no more right clicks)
    % time of next left click (set to infinity if no more left clicks)
    % stimulusLength
% update vector index of click train, if used
% evolve ODE exactly up to t_new
% add delta jump if t_new is the time of a click

while time<stimulusLength
    [t_new,ii] = min([lnxt, rnxt, stimulusLength]);
    if ii == 1
        jump = -kappa;
        ltidx = ltidx + 1;
        if ltidx > lmax
            lnxt = inf;
        else
            lnxt = lTrain(ltidx);
        end
    elseif ii == 2
        jump = kappa;
        rtidx = rtidx + 1;
        if rtidx > rmax
            rnxt = inf;
        else
            rnxt = rTrain(rtidx);
        end
    else
        jump=nan;
    end         
    
    %evolve ODE exactly between time and t_new according to equation
    % y(time)=a
    % y(t_new)=2*acoth(coth(a/2)*exp(2*h*(t_new-time)))
    evidence=2*acoth(coth(evidence/2)*exp(2*h*(t_new-time)));
    if not(isnan(jump))
        evidence = evidence + jump;
    end
    time = t_new;
end
ss=sign(evidence);
end