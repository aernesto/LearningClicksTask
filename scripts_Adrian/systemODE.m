function [ss,tt,uu]=systemODE(lTrain, rTrain, rateLow, rateHigh, T, gammax,...
posttimes, priorState, alpha, beta, dt, cptimes)
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
%   dt - timestep to use for Euler method
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
fileID=fopen('sysODElog.txt','w');
while time<T
    jump=0;
    t_new=time + dt;
    inttime=floor(t_new);
    % print to log file
    %if and(abs(t_new-inttime)<2*dt, inttime == 27)
     %   fprintf(fileID,'t_new = %.3f \n',t_new);
    %end
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
    
    % following if statement is probably redundant
    if not(presclick)
        jump=0;
    end         
    presclick=false;
    
    % evolve log posterior vectors
        % scalar boundary case for gamma=0
    %%%%%%%%%%%%%%%%
    % BUG NEXT LINE
    % returns NaN values
    %
    % following is an extract from log file:
    %%%
    % t_new = 25.9747 
    % length(yp_new) = 100 
    % isnan(yp_new) = 0 
    % t_new = 25.9748 
    % length(yp_new) = 100 
    % isnan(yp_new) = 1 
    %%%
    %%%%%%%%%%%%%%%%
    yp_new_gamma0=yp_old(1)-alpha*dt/(time+beta);
    ym_new_gamma0=ym_old(1)-alpha*dt/(time+beta);
    %if and(abs(t_new-inttime)<2*dt, inttime == 27)
     %   fprintf(fileID,'yp_old(1) = %.5f \n', yp_old(1));
      %  fprintf(fileID,'isnan(yp_old(1)) = %d \n', isnan(yp_old(1)));
        %fprintf(fileID,'alpha-1 = %.3f \n', alpha-1);
       % fprintf(fileID,'exp(-yp_old(1)) = %.3f \n', exp(-yp_old(1)));
        %fprintf(fileID,'dt/(time+beta) = %.3f \n', dt/(time+beta));
       % fprintf(fileID,'yp_new_gamma0 = %.3f \n',yp_new_gamma0);
       % fprintf(fileID,'ym_new_gamma0 = %.3f \n',ym_new_gamma0);
    %end
        % rest of vectors
    yp_prime=-(alpha+gammaValues(2:end))';
    ym_prime=yp_prime;
    
    yp_prime=yp_prime+(alpha-1+gammaValues(2:end)').*exp(ym_old(1:end-1)-yp_old(2:end));
    %if and(abs(t_new-inttime)<2*dt, inttime == 27)
      %      fmt = [repmat('%4d ', 1, size(yp_prime,2)-1), '%4d\n'];
     %       fprintf(fileID,'NaN yp_prime = \n');
            %fprintf(fileID,fmt,isnan(yp_prime).');
    %end
    
    yp_prime=yp_prime / (time + beta);
    
    ym_prime=ym_prime+(alpha-1+gammaValues(2:end)').*exp(yp_old(1:end-1)-ym_old(2:end));
    ym_prime=ym_prime / (time + beta);
    
        % concatenate
    yp_new = [yp_new_gamma0;dt*yp_prime + yp_old(2:end)];
    ym_new = [ym_new_gamma0;dt*ym_prime + ym_old(2:end)];
    
    % add jump 
    yp_new = yp_new + jump;
    ym_new = ym_new + jump;
    
    % if report time hit, normalize and output posterior variance
    if t_new > nxtposttime
        %normalization constant
        
        
        % on trial 2: ERROR
        %t_new = 27.62340 
        %K = -739.23607 
        %t_new = 27.62350 
        %K = -Inf
        %
        %
        K = log(sum(exp(yp_new)+exp(ym_new)));
        %true posterior
        xp=yp_new-K;
        xm=ym_new-K;
    
        %posterior variance over h
        v1=(gammaValues'+alpha)/(t_new+beta);
        v2=(gammaValues'+alpha+1)/(t_new+beta);
        AVG=sum((exp(xp)+exp(xm)).*v1);
        post_mean_h(idnxtposttime)=AVG;
        post_var_h(idnxtposttime)=sum((exp(xp)+exp(xm)).*v1.*v2)-...
            sum((exp(xp)+exp(xm)).*v1)^2;
        
        %report lower bound on posterior variance
        %(alpha+n)/(beta+t)^2
        ncp = sum(cptimes<t_new); % number of true change points by report time
        lbvar(idnxtposttime)=(alpha+ncp)/(beta+t_new)^2;
        
        if abs(inttime-27) < 3
            fprintf(fileID,'----------------------------\n');
            fprintf(fileID,'report time report \n');
            fprintf(fileID,'t_new = %.5f \n',t_new);
            fprintf(fileID,'mean = %.5f \n',AVG);
            fprintf(fileID,'----------------------------\n');
            %fmt = [repmat('%4d ', 1, size(yp_new,2)-1), '%4d\n'];
            %fprintf(fileID,'NaN yp_new = \n');
            %fprintf(fileID,fmt, isnan(yp_new).');
            
            %disp(isnan(yp_new'));
            %fprintf(fileID,'NaN ym_new \n');
            %disp(isnan(ym_new));
            
            %fprintf(fileID,'nxtposttime = %.3f \n', nxtposttime);
            %fprintf(fileID,'K = %.3f \n',K);
            
            %fprintf(fileID,'NaN v1 \n');
            %disp(isnan(v1'));
            %fprintf(fileID,'NaN v2 \n');
            %disp(isnan(v2));
            %fprintf(fileID,'NaN mean \n');
            %disp(isnan(post_mean_h(idnxtposttime)));
            %fprintf(fileID,'NaN varexp \n');
            %disp(isnan((exp(xp)+exp(xm)).*v1.*v2)');
        end
        
        
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
    
    if (27 < t_new) && (t_new < 27.7)
        fprintf(fileID,'t_new = %.5f \n',t_new);
        fprintf(fileID,'K = %.5f \n',log(sum(exp(yp_new)+exp(ym_new))));
        %fprintf(fileID,'isnan(ym_new) = %d \n',sum(ym_new==Inf));
    end
    
    time = t_new;
end
fclose(fileID);
tt=post_mean_h;
ss=post_var_h;
uu=lbvar;
end