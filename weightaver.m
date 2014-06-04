function f_est = weightaver( freqs,segSNR )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

f_est=sum(freqs.*segSNR)/sum(segSNR);
end

