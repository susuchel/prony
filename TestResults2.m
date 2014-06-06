%f_est=[];
% analysis of results

f0;
flim=[10000 15000 20000 28500 30000 34000 40000 50000]; % length = length(f0)+1
C=cell(length(f0),4);
for k=1:length(f0) % k stands for component to be estimated
    f_est=[];f_mean=[]; F=[]; F_s=[]; df_est=[]; df_mean=[]; 
    F_s = choose_freq(f,segSNR,flim(k),flim(k+1));
    %[f_est,F,F_s]=weightaver(f,segSNR,flim(k),flim(k+1));
    f_est=weightaver(F_s(:,1),F_s(:,2),7);
    for l=1:length(F_s(:,1)) 
        f_mean(l,1)=mean(F_s(1:l,1)); 
    end
    
    df_est=abs(f_est-f0(k))*100/f0(k);
    df_mean=abs(f_mean-f0(k))*100/f0(k);
    
    figure(20+k)
    subplot(2,1,1); plot([f_est,f_mean,f0(k)*ones(length(f_est),1)])
    subplot(2,1,2); plot([df_est,df_mean])
    
    C{k,1}=f_est; C{k,2}=df_est; C{k,3}=f_mean; C{k,4}=df_mean;
    
end