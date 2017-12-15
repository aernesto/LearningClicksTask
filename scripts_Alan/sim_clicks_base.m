% core code to simulate click task 
% Alan Veliz-Cuba 12/11/17
clear;
h=1; %hazard rate
rateLow=2;
rateHigh=38;
kappa=log(rateHigh/rateLow); 
stimulusLength=0.5;
numTrials=10000;

%the following may have to be modified in Matlab
seed=1; randn('state',seed); rand('state',seed)
dt=0.001; 
%interval from 0 to stimulusLength
INT=0:dt:stimulusLength;
%number of iterations
num_it=length(INT);
%initialize y_t and environment
y=zeros(numTrials,num_it);
E=zeros(numTrials,num_it);
%create clicks from the high and low rate
%the assignment to left and right is done later
highClicks=rand(numTrials,num_it)<rateHigh*dt;
lowClicks=rand(numTrials,num_it)<rateLow*dt;
%initial environment
E(:,1)=1;
tic
%define -dt*2*h to speed up sims
minus2hdt=-dt*2*h;
for i=1:num_it-1
    %update environment. Note: E\in {-1,1}
    E(:,i+1)=E(:,i).*(-1).^(rand(numTrials,1)<h*dt);
    %update y using simple Euler's method
    y(:,i+1)=y(:,i) + minus2hdt*sinh(y(:,i));
    %update y using clicks. the mult by E(:,i) assigns left and right
    y(:,i+1)=y(:,i+1)+kappa*(highClicks(:,i)-lowClicks(:,i)).*E(:,i);%.*(-1).^(rand(numTrials,1)<q);
end

%get the sign of y_t
signY=sign(y);
%the following forces the prediction to be random if y=0
prediction=signY.*abs(signY)+(1-abs(signY)).*(-1).^(rand(numTrials,num_it)<.5);
%compute the average accuracy across trials
correct=mean(prediction==E);
%plotting
plot(INT,correct,'LineWidth',2); 
xlabel('interrogation time')
ylabel('percentage correct')
title(['h=' num2str(h) ', lambda_L=' num2str(rateLow) ', lambda_H=' num2str(rateHigh)]);
