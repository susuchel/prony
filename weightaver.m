function f_est = weightaver( freqs,segSNR,f2,f1 )
%weightaver Takes weighted average of f estimate in MProny
%   Detailed explanation goes here

% TODO: Add deleting zero frequencies, complete this function
% TODO: Create new function to estimate error for ployharmnic  signal 
kk=1;
freqs1=0;
segSNR1=0;
for k=1:length(freqs)
    %if freqs(k)~=0
    if ((freqs(k)<f2)&&(freqs(k)>f1))
        freqs1(kk)=freqs(k);
        segSNR1(kk)=segSNR(k);
        kk=kk+1;
    end
end
f_est=sum(freqs1.*segSNR1)/sum(segSNR1);
end

