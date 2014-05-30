function [sp,ffs] = spfft (fs,signal,fig_mode)
% spfft (fs,signal,fig_mode) FFT with required parameters
%   fs - sampling frequency 
%   signal - signal to analyze, must be a column, if signal is a matrix, 
%   then fft will be applied to each coloumn
%   fig_mode - variable that tells whether to plot or not to plot 
%   if fig_mode=0 or =false, then it do not plot a graph

qwe=signal;
[a,b]=size(qwe); % a is number of rows, b is number of columns
% each column should be an implementation of the signal
% a is number of samples in an implementation
meanqwe=mean(qwe);
for k=1:b
    qwe(:,k)=qwe(:,k)-meanqwe(1,k);
end
    
Nfft = a; % must be even

y=fft(qwe,Nfft); 
Pyy = y.* conj(y)/Nfft;
ff = (fs*(0:(Nfft-1))/Nfft)';

sp=Pyy(1:Nfft/2,:); 
ffs=ff(1:Nfft/2)/1000;


if fig_mode ~=false
    spnorm=zeros(size(sp));
    spmax=max(sp);
    for k=1:size(sp,2) 
        spnorm(:,k)=sp(:,k)/max(1,k);
    end
    figure(1);
    plot(ffs/1000,sp/max(sp));
end
end