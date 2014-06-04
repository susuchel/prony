function SNR = calc_SNR2 (signal,sig)
%calc_SNR2 Calculates signal-to-noise ratio
%   Uses estimate of std of noise
%
%   Inputs:
%   signal  - signal to analyze, may be a matrix, where each column is an 
%   observation 
%   sig     - std of noise, or its unbiased estimate in case of a real signal
%   SNR = mean(10*log10(sum(signal.*signal)/(sig*sig*size(signal,1))-1));
if (sig<=0)
    error('calc_SNR2 (signal,sig): ''sig'' must be positive')
end
SNR = 10*log10(sum(signal.*signal)/(sig*sig*size(signal,1))-1);
SNR(real(SNR)<0)=0;
end

 