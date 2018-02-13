function ss=returnPostH(lTrain, rTrain, rateLow, rateHigh, T, gammax,...
posttimes, priorState, alpha, beta, dt, cptimes)
% DESCRIPTION:
% This function evolves the system of jump ODEs for the unknown hazard 
% rate, over a given stimulus train, and returns the values of the
% posterior marginal over CP count at the prescribed assessment times.
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
%   dt - timestep to use for Euler method
%
% RETURNS:
%   ss - stochastic matrix (cols sum to 1) of dimensions 
%        2*gammax-by-posttimes containing the 
%        values of the joint posterior over state and CP count

ss=zeros(2*gammax,length(posttimes));
coliter=1;
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
%priorGamma=((alpha.^gammaValues)*exp(-alpha))./factorial(gammaValues);
massOn0=.99;
priorGamma=[massOn0,ones(1,gammax-1)*(1-massOn0)/gammax];
yp_old=log(priorState(1)*priorGamma)';
ym_old=log(priorState(2)*priorGamma)';
time=0;

% set delta jump sizes
kappa = log(rateHigh/rateLow);

post_var_h=zeros(size(posttimes));
post_mean_h=post_var_h;
lbvar=post_var_h;%vector storing theoretical lower bound on variance
nposttimes = length(posttimes);
if nposttimes > 0
    idnxtposttime = 1;
    nxtposttime=posttimes(1);
else
    nxtposttime=inf;
end

%Forward Euler
presclick=false; % if true, means at least one click fell in current time bin
%fileID=fopen('sysODElog.txt','w');
while time<T
    jump_right=0;
    jump_left=0;
    t_new=time + dt;
    inttime=floor(t_new);
    if t_new > lnxt
        presclick=true;
        jump_left = kappa;
        ltidx = ltidx + 1;
        if ltidx > lmax
            lnxt = inf;
        else
            lnxt = lTrain(ltidx);
        end
    end
    if t_new > rnxt
        presclick=true;
        jump_right = kappa;
        rtidx = rtidx + 1;
        if rtidx > rmax
            rnxt = inf;
        else
            rnxt = rTrain(rtidx);
        end
    end
    
    % following if statement is probably redundant
    if not(presclick)
        jump=0;
    end         
    presclick=false;
    
    % evolve log posterior vectors
        % scalar boundary case for gamma=0
    yp_new_gamma0=yp_old(1)-alpha*dt/(time+beta);
    ym_new_gamma0=ym_old(1)-alpha*dt/(time+beta);
        % rest of vectors
    yp_prime=-(alpha+gammaValues(2:end))';
    ym_prime=yp_prime;
    
    yp_prime=yp_prime+(alpha-1+gammaValues(2:end)').*exp(ym_old(1:end-1)-yp_old(2:end));
    
    yp_prime=yp_prime / (time + beta);
    
    ym_prime=ym_prime+(alpha-1+gammaValues(2:end)').*exp(yp_old(1:end-1)-ym_old(2:end));
    ym_prime=ym_prime / (time + beta);
    
        % concatenate
    yp_new = [yp_new_gamma0;dt*yp_prime + yp_old(2:end)];
    ym_new = [ym_new_gamma0;dt*ym_prime + ym_old(2:end)];
    
    % add jump 
    yp_new = yp_new + jump_right;
    ym_new = ym_new + jump_left;
    
    % if report time hit, normalize and output posterior variance
    if t_new > nxtposttime
        %normalization constant
        K = log(sum(exp(yp_new)+exp(ym_new)));
        %true posterior
        xp=yp_new-K;
        xm=ym_new-K;
        vector=[exp(xp);exp(xm)];
        %disp(t_new)
        %disp(sum(marginalOverGamma))
        ss(:,coliter)=vector;
        coliter=coliter+1;
               
        %update index of next reporting time
        if idnxtposttime == nposttimes
            nxtposttime = inf;
        else
            idnxtposttime = idnxtposttime + 1;
            nxtposttime = posttimes(idnxtposttime);
        end
        yp_old = xp;
        ym_old = xm;
    else
        yp_old = yp_new;
        ym_old = ym_new;
    end
    
    % reinitialize for next iteration    
    time = t_new;
end
end