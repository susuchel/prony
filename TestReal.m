clearvars -except T_RES
tic;
% load a file with data
filename=input('Please define the filename using '' '':\n');
load(filename, 'signal', 'fs', 'sig', 'sp', 'ffs')

% main inputs are in file: 
% signal    - signal to analyze (rows - samples, columns - observations)
% fs        - sampling frequency
% sig       - std of noise of each observation
% sp        - spectrum of the signal
% ffs       - array of corresponding frequecies to sp


% plot signal and its FFT spectrum
spnorm=zeros(size(sp)); spmax=max(sp);
for k=1:size(sp,2) 
    spnorm(:,k)=sp(:,k)/spmax(1,k); 
end

figure(1); plot(signal)
figure(2); plot(ffs/1000, spnorm)

% define flim
% Note that flim for real signals is in another form
flim=input('Please define flim in Hz in the form: [lim1_start,lim1_end;...;limn_start,limn_end] \n');

% fft estimates of frequencies, chosen frequencies to control
f0_fft=choose_fft(sp,ffs,flim); 

% Fs_opt
fopt=find_fopt(f0_fft); % you should check whether fopt>=2*max(f0)
                            % otherwise consider filtering
                            % Note that there may be frequencies greater than
                            % max(flim), so be careful when using resampling 

str1='Chosen FFT frequencies are ';
str2=sprintf('%g ',f0_fft);
str3=sprintf('\nOptimal frequency is %g\nChange fs? If yes, enter 1, if no, enter 0\n',fopt);
str=[str1,str2,str3];
change_fs_mode=input(str);
if change_fs_mode~=0
    [signal, fs]=change_fs(signal,fs,fopt);                      
end


% MProny inputs

p=input('Specify p = ');
N_samp=input('Specify length of a segment N_samp = ');
eps=input('Specify eps = ');

fig_mode=0; 
hist_mode=1; 
savefilename=0;

svd_mode=input('Use svd? If yes, enter 1, if no, enter 0\n'); 
if svd_mode~=0
    q=input('Specify q = ');
end

% MProny
obs_lim=input('Please define obs_lim in the form n1:n2\n');
[freqs,alphas,segSNR,segSV]=MProny(signal(:,obs_lim),fs,sig,p,q,N_samp,eps,svd_mode,fig_mode,hist_mode,savefilename);


hlim=zeros(1,2*size(flim,1)); % for histogram in our limits 
for k=1:size(flim,1) 
    hlim(2*k-1)=flim(k,1); 
    hlim(2*k)=flim(k,2);
end
figure(3); hist(freqs,hlim) % Note that all limits in flim must be different

% ? plot(segSNR) - to choose right SNRlevel
% ? plot(segSV) - to choose p or q

% ? ask to change flim
SNRlevel=input('Specify SNRlevel = ');

T_res=table();
for k=1:size(flim,1)
    Fsorted = choose_freq(freqs,segSNR,flim(k,1),flim(k,2));
    [f_est,f_mean]=weightaver(Fsorted(:,1),Fsorted(:,2),SNRlevel);
    df_est=(abs(f_est-f0_fft(k))*100/f0_fft(k))'; 
    df_mean=(abs(f_mean-f0_fft(k))*100/f0_fft(k))';
    % TODO: consider using cell arrays with tables as cells
    T_temp=table(flim(k,:), f0_fft(k), f_est(end), df_est(end),f_mean(end), df_mean(end));
    T_res=[T_res;T_temp];
end 

toc;