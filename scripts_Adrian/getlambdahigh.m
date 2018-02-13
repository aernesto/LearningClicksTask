function rates=getlambdahigh(lambdalow, snr, upperroot)
% DESCR: this function returns the values of lambdahigh required to reach
% the prescribed values of SNR, for a given lambdalow value. The formula
%  snr=(lambdahigh-lambdalow)./sqrt(lambdahigh+lambdalow);
% is inverted, taking the root of the quadratic specified by upperroot.
% ARGS:
%   lambdalow= positive scalar in Hz (click rate)
%   snr= row vector of desired SNR values (positive values)
%   upperrot = TRUE if result should be upper root; FALSE otherwise
% RETURNS:
%    row vector of lambdahigh values for each snr value
if upperroot
    rates=0.5*(2*lambdalow+snr.^2+sqrt(snr.^4+4*lambdalow*snr.^2+4*lambdalow*snr.^2));
else
    rates=0.5*(2*lambdalow+snr.^2-sqrt(snr.^4+4*lambdalow*snr.^2+4*lambdalow*snr.^2));
end
end