function [f_est, f_mean]= weightaver(freqs,segSNR,SNRlevel)
%weightaver Takes weighted average of f estimate in MProny
%   TODO: Detailed explanation goes here
%   freqs should be one column
%   output is frequency estimate depending on number of segments to take
%   weighted average
%   Note that in this function data with SNR=0 or SNR<SNRlevel does not 
%   change result, because fi*0=0

f_est=[];
f_mean=[];
if (SNRlevel~=0)
    segSNR(segSNR<SNRlevel)=0;
end
for k=1:length(segSNR)
    f_est(k,1)=sum(freqs(1:k,1).*segSNR(1:k,1))/sum(segSNR(1:k,1));
    good=segSNR(1:k,1)>=SNRlevel;
    f_mean(k,1)=mean(freqs(good,1));
end
end

