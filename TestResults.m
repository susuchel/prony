clear all
clc
tic;
count=0;
%profile on
% changing inputs
SNR_val=5;%[15, 10, 5];%[40,20:-5:5];
N_obs_max=10;%3;
N_obs_val=1:N_obs_max;
eps_val=0.05;
p_val=64;%[32,64];%[4,8,16,32,64];
N_samp_val=200; %[200,300];%[128, 200, 300, 400, 500];
SNRlevel_val=5;%[0 5 10];


% fixed inputs 
fs=500000; 
f0=[13853, 17690, 27690, 29600, 32600, 35421, 44310]; 
alfa=[610, 2000, 1060, 1900, 1480, 870, 3200]; 
N_mod=2500; 
flim=[10000, 15000, 20000, 28500, 31500, 33000, 40000, 48000];
fig_mode_fft=0; fig_mode_MP=0; hist_mode=0; filename=0;

svd_mode=1; q=2*length(f0);

change_fs_mode=1;

T=table(); % for results
%T2=table();


for SNR=SNR_val
    
    % model generation
    [m,s]=model1(fs,f0,alfa,N_mod,SNR,N_obs_max);

    if change_fs_mode~=0
        fopt=find_fopt(f0);
        [m, fs_new]=change_fs(m,fs,fopt);
    else
        fs_new=fs;
    end
    
    %for fs=[fs,fs_new]
        
        for N_obs=N_obs_val

            % spfft
            [sp,ffs]=spfft(fs_new,m(:,1:N_obs),fig_mode_fft);
            f0_fft = choose_fft(sp,ffs,flim);
            df_fft=abs(f0_fft-f0)*100./f0; 

            for eps=eps_val
                for p=p_val
                    for N_samp=N_samp_val
                        % MProny
                        [f,a,segSNR,segSV]=MProny(m(:,1:N_obs),fs_new,s,p,q,N_samp,eps,svd_mode,fig_mode_MP,hist_mode,filename);

                        for k=1:length(f0)

                            F_s = choose_freq(f,segSNR,flim(k),flim(k+1));

                            for SNRlevel=SNRlevel_val

                                [f_est,f_mean]=weightaver(F_s(:,1),F_s(:,2),SNRlevel);

                                df_est=(abs(f_est-f0(k))*100/f0(k))'; 
                                df_mean=(abs(f_mean-f0(k))*100/f0(k))';

                                % result 
                                count=count+1;
                                
                                C{count,1}=f_est;
                                C{count,2}=df_est';
                                C{count,3}=f_mean;
                                C{count,4}=df_mean';
                                
                                maxdf_est=max(df_est); 
                                maxdf_mean=max(df_mean);
                                
                                enddf_est=df_est(end);
                                enddf_mean=df_mean(end);

                                T2=table(fs_new,f0,alfa,N_mod,flim,svd_mode,q,SNR,N_obs,df_fft,eps,p,N_samp,[flim(k),flim(k+1)],f0(k),SNRlevel,maxdf_est,maxdf_mean,enddf_est,enddf_mean);
                                T=[T;T2];
                                
                                
                            end
                        end
                    end
                end
            end
        end
    %end
end

T.Properties.VariableNames={'fs','f0','alfa','N_mod','flim','svd_mode','q','SNR','N_obs','df_fft','eps','p','N_samp','flim_local','f0_local','SNRlevel','maxdf_est','maxdf_mean','enddf_est','enddf_mean'};
toc;
count
%profile viewer