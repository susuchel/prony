% TODO: add explanations

flim=[13000,15000;
      16000,17000;
      17000,18500;
      26000,27000;
      29000,30000;
      32000,33000;
      36500,38000;
      44000,45000;
      45200,46000;
      57000,59000];
      %63000,65000];
  
change_fs_mode=1;
fopt=135000;

C_res=[];
f0_test=16542;
m=model1(500000,f0_test,610,2500,0,500); % test model to add to signal
% take signal
for k_file=1:8
    filename=sprintf('Yliq%d.mat',k_file);
    load(filename, 'signal', 'fs', 'sig')
    
    for kk=1:size(signal,2) 
        signal(:,kk)=signal(:,kk)+m; 
    end
    C_res{k_file,1}=filename;
        
    for k=1:size(flim,1)
        D_res{k,1}=flim(k,:);
        D_res{k,2}=table();
    end
    
    if change_fs_mode~=0
        %fopt=find_fopt(f0);
        [signal, fs]=change_fs(signal,fs,fopt);
    end
          
    for N_obs=10
              
        [sp,ffs]=spfft(fs,signal(:,1:N_obs),0);
        f0_fft=choose_fft(sp,ffs,flim);
        
      for len=625 % to decrease number of segments
        
        for eps=[0.1, 0.05]
            fig_mode=0; 
            hist_mode=0; 
            savefilename=0;

            for p=[80,120,144]
                N_samp=300;
                
                for svd_mode=1
                    q=24;
                    [freqs,alphas,segSNR,segSV]=MProny(signal(1:len,1:N_obs),fs,sig,p,q,N_samp,eps,svd_mode,fig_mode,hist_mode,savefilename);
                    
                    for SNRlevel=5
                        
                        for k=1:size(flim,1)
                            Fsorted = choose_freq(freqs,segSNR,flim(k,1),flim(k,2));
                            [f_est,f_mean]=weightaver(Fsorted(:,1),Fsorted(:,2),SNRlevel);
                            df_est=(abs(f_est-f0_fft(k))*100/f0_fft(k))'; 
                            df_mean=(abs(f_mean-f0_fft(k))*100/f0_fft(k))';
                            if flim(k,1)==16000
                                df_est=(abs(f_est-f0_test)*100/f0_test)'; 
                                df_mean=(abs(f_mean-f0_test)*100/f0_test)';
                            end
    
                            T_temp=table(fs,f0_fft(k),f_est(end), df_est(end),f_mean(end),df_mean(end),SNRlevel,svd_mode, q, p, N_samp,eps, N_obs,len);
                            D_res{k,2}=[D_res{k,2};T_temp];
                        end
                    end
                end
            end
        end
      end
    end
    
    for k=1:size(flim,1)
        D_res{k,2}.Properties.VariableNames={'fs','f0_fft','f_est','df_est','f_mean','df_mean','SNRlevel','svd_mode', 'q', 'p', 'N_samp','eps', 'N_obs','len'};
    end
        
    C_res{k_file,2}=D_res;
end

RESopt3=C_res;
save('Data\RESopt.mat','RESopt3','-append')

