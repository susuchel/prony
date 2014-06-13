%f_est=[];
% analysis of results

%flim=[10000 15000 20000 25000 30000 35000 40000 45000 60000 70000]; % length = length(f0)+1

for k=1:(length(flim)-1) % k stands for component to be estimated
    f_est=[];f_mean=[]; F_s=[]; df_est=[]; df_mean=[]; 
    F_s = choose_freq(f144svd(:,:,1),segSNR(:,1),flim(k),flim(k+1));
    %ff0=ffs(sp(:,1)==max(sp(flim(k)/200:flim(k+1)/200,1)),1)
    %f_est=weightaver(F_s,segSNR(:,1),flim(k),flim(k+1));
    f_est=weightaver(F_s(:,1),F_s(:,2),0);
    for l=1:length(F_s(:,1)) 
        f_mean(l,1)=mean(F_s(1:l,1)); 
    end
    
%     df_est=abs(f_est-f0(k))*100/f0(k);
%     df_mean=abs(f_mean-f0(k))*100/f0(k);
    
    figure(k+20)
%     subplot(2,1,1); plot([f_est,f_mean,f0(k)*ones(length(f_est),1)])
%     subplot(2,1,2); plot([df_est,df_mean])
    plot([f_est,f_mean,f0_fft(k)*ones(length(f_est),1), (f0_fft(k)+200)*ones(length(f_est),1), (f0_fft(k)-200)*ones(length(f_est),1)])
    
    C{k,1}=f_est; C{k,2}=f_mean; C{k,3}=F_s;
    
end