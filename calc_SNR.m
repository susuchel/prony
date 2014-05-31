function SNR = calc_SNR (signal)
%calc_SNR Calculates signal-to-noise ratio
%   Try to calculate SNR based on fft. 
%   Does not work properly.
%   Estimates are acceptable only about 6 dB

Nfft = size(signal,1); % must be even
y=fft(signal,Nfft); 
Pyy = y.* conj(y)/Nfft;
a=(Pyy>mean(mean(Pyy)));
b=Pyy.*a;
n=Pyy-b;
SNR=mean(10*log10(sum(b)./sum(n)));
end

