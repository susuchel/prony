% clear all
% clc
% tic;
%count=0;
%profile on
% changing inputs
SNR_val=[10, 5];%[40,20:-5:5];
N_obs_max=10;%3;
N_obs_val=10;%1:N_obs_max;
eps_val=[0.1 0.05 0.01];
p_val=[64,80,100,120];%[32,64];%[4,8,16,32,64];
N_samp_val=[250,300,400]; %[200,300];%[128, 200, 300, 400, 500];
SNRlevel_val=5;%[0 5 10];


% fixed inputs 
fs=500000; 
f0=[13853, 17690, 27690, 29600, 32600, 35421, 44310]; 
alfa=[610, 2000, 1060, 1900, 1480, 870, 3200]; 
N_mod=2500; 
flim=[13000, 15000; 
      17000,18500; 
      25000, 28500; 
      29000, 31500; 
      32000, 33000; 
      34000,37000; 
      43000, 45000];
fig_mode_fft=0; fig_mode_MP=0; hist_mode=0; filename=0;

svd_mode=1; q=2*length(f0);

change_fs_mode=1;

%T=table(); % for results
%T2=table();

D_res=[];
for k=1:size(flim,1)
    D_res{k,1}=flim(k,:);
    D_res{k,2}=table();
end


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
                        len=500;
                        % MProny
                        [f,a,segSNR,segSV]=MProny(m(1:len,1:N_obs),fs_new,s,p,q,N_samp,eps,svd_mode,fig_mode_MP,hist_mode,filename);

                        for k=1:size(flim,1)

                            F_s = choose_freq(f,segSNR,flim(k,1),flim(k,2));

                            for SNRlevel=SNRlevel_val

                                [f_est,f_mean]=weightaver(F_s(:,1),F_s(:,2),SNRlevel);

                                df_est=(abs(f_est-f0(k))*100/f0(k))'; 
                                df_mean=(abs(f_mean-f0(k))*100/f0(k))';

%                                 maxdf_est=max(df_est); 
%                                 maxdf_mean=max(df_mean);
%                                 
%                                 enddf_est=df_est(end);
%                                 enddf_mean=df_mean(end);

                                T_temp=table(fs_new,SNR, alfa(k),f0(k),f0_fft(k), df_fft(k),f_est(end), df_est(end),f_mean(end),df_mean(end),SNRlevel,svd_mode, q, p, N_samp,eps, N_obs,len);
                                D_res{k,2}=[D_res{k,2};T_temp];
%                                 T2=table(fs_new,f0,alfa,N_mod,svd_mode,q,SNR,N_obs,df_fft,eps,p,N_samp,flim(k,:),f0(k),SNRlevel,maxdf_est,maxdf_mean,enddf_est,enddf_mean);
%                                 T=[T;T_temp];
                                
                                
                            end
                        end
                    end
                end
            end
        end
    %end
end

for k=1:size(flim,1)
    D_res{k,2}.Properties.VariableNames={'fs_new','SNR','alfa', 'f0','f0_fft', 'df_fft','f_est', 'df_est','f_mean','df_mean','SNRlevel','svd_mode', 'q', 'p', 'N_samp','eps', 'N_obs','len'};
end

%profile viewer
RES_MOD2=D_res;
save('RES_MOD.mat','RES_MOD2','-append')