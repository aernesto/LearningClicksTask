function ss=systemODE(lTrain, rTrain, rateLow, rateHigh, T, gammax,...
posttimes, priorState, alpha, beta)
% DESCRIPTION:
% This function evolves the system of jump ODEs for the unknown hazard 
% rate, over a given stimulus train, and returns the values of the
% posterior variance of h at the prescribed assessment times.
%
% ARGUMENTS:
%   lTrain - column vector of left train click times in sec
%   rTrain - column vector of right train click times in sec
%   rateLow - low click rate in Hz
%   rateHigh - high click rate in Hz
%   T - trial duration
%   gammax - maximum number of change point counts allowed.
%               note that gammax=50 means that -1<gamma<50
%   posttimes - column vector of times at which the posterior variance
%               should be computed
%   priorState - 2x1 vector containing the prior probabilities over the 
%                   environmental states H+ and H- respectively
%   alpha - hyperparameter of Gamma dist over h
%   beta - hyperparameter of Gamma dist over h
%
% RETURNS:
%   ss - column vector with same dimensions as posttimes containing the 
%           values of the posterior variance over hazard rate

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

% initialize posterior vectors and time
gammaValues=0:gammax-1;
% Poisson prior over change point counts
priorGamma=((alpha.^gammaValues)*exp(-alpha))./factorial(gammaValues);
yp_old=log(priorState(1)*priorGamma)';
ym_old=log(priorState(2)*priorGamma)';
time=0;

% set delta jump sizes
kappa = log(rateHigh/rateLow);

%Forward Euler
presclick=false; % if true, means at least one click fell in current time bin
while time<stimulusLength
    jump=0;
    t_new=time + dt;
    if t_new > lnxt
        presclick=true;
        jump = -kappa;
        ltidx = ltidx + 1;
        if ltidx > lmax
            lnxt = inf;
        else
            lnxt = lTrain(ltidx);
        end
    end
    if t_new > rnxt
        presclick=true;
        jump = jump+kappa;
        rtidx = rtidx + 1;
        if rtidx > rmax
            rnxt = inf;
        else
            rnxt = rTrain(rtidx);
        end
    end
    
    %[t_new,ii] = min([lnxt, rnxt, stimulusLength]);
    if not(presclick)
        jump=0;
    end         
    presclick=false;
    
    % evolve posterior vectors
        % scalar boundary case for gamma=0
    yp_new_gamma0=yp_old(1)+((alpha-1)*exp(-yp_old(1))-alpha)*dt/(time+beta);
    ym_new_gamma0=ym_old(1)+((alpha-1)*exp(-ym_old(1))-alpha)*dt/(time+beta);
    
        % rest of vectors
        
        %%%%%%%%%%%%%%%%%%%%%
        %%%%%% TAKE UP HERE %%%%%%%%%%%%%%%
    yp_new=yp_old(2:end)+((alpha-1)*exp(-yp_old(1))-alpha)*dt/(time+beta);
    
        % concatenate
    yp_new = [yp_new_gamma0;yp_new];
    ym_new = [ym_new_gamma0;ym_new];
    % add jump 
    yp_new = yp_new + jump;
    ym_new = ym_new + jump;
        
    time = t_new;
end
ss=sign(evidence);
end