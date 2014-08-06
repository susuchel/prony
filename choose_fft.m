function f0_fft = choose_fft(sp,ffs,flim)
%CHOOSE_FFT Determines frequency estimates by FFT within intervals
%   choose_fft(spfft,ffs,flim) returns an array of frequency estimates 
%   within intervals managed by flims
%
%   Inputs:
%   sp      - power spectrum from spfft function
%   ffs     - corresponding frequencies from spfft function
%   flim    - limits of intervals to find maximum. 
%           FIXME: add new description of flim structure
%           OLD:It has following structure:
%           [flim(1),flim(2),flim(3),...,flim(n)]. In this case function will
%           return (n-1) frequencies, corresponding local maximum within
%           the intervals with limits flim(k) and flim(k+1). In flim one
%           should specify frequencies in Hz
%

%TODO:  function will return f within the interval even if relative value of
%       maximum within the interval with respect to values within other intervals 
%       is very small. Need to be compared with some parameter, or function 
%       may return not only frequencies but corresponding values of sp 
% FIXME: change description, because now it take average over observations

df=ffs(2)-ffs(1);
[~,N_obs]=size(sp);
f0_fft_obs=zeros(size(flim,1),N_obs);
for i_obs=1:N_obs
    for k=1:size(flim,1)
        fftemp=ffs(uint16(flim(k,1)/df+1):uint16(flim(k,2)/df+1),1);
        f0_fft_obs(k,i_obs)=fftemp(sp(uint16(flim(k,1)/df+1):uint16(flim(k,2)/df+1),i_obs)==max(sp(uint16(flim(k,1)/df+1):uint16(flim(k,2)/df+1),i_obs)),1);
    end
end
if N_obs~=1
    f0_fft=mean(f0_fft_obs'); % take mean over obseravtions
else
    f0_fft=f0_fft_obs';
end

