function SNR = calc_SNR2 (signal,sig)
%calc_SNR2 Calculates signal-to-noise ratio
%   Uses estimate of std of noise
%
%   Inputs:
%   signal  - signal to analyze
%   sig     - std of noise, or its unbiased estimate in case of a real signal
SNR = mean(10*log10(sum(signal.*signal)/(sig*sig*size(signal,1))-1));
end

 