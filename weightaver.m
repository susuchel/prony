function f_est = weightaver( freqs,segSNR )
%weightaver Takes weighted average of f estimate in MProny
%   Detailed explanation goes here

f_est=sum(freqs.*segSNR)/sum(segSNR);
end

