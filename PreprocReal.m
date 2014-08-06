% preproccessing of a real signal
% aim is to define parameters of signal

%clearvars -except Yout
clear all

% choose signal and define fs
% be sure that you add a directory with file to path 

filename=input('Please define the filename using '' '':\n');
load(filename);
signal=Yout'; % each column should be an observation
fs=500000;

fprintf('The signal %s has %d obseravations, each of %d samples\n', filename, size(signal,2), size(signal,1)); 

% first you should define informative part of the signal, 
% do that using picture
% N_obs_plot=input('How many obs to plot?\n');
N_obs_plot=5000;
if N_obs_plot~=0 
    figure(1); plot(signal(:,1:N_obs_plot))
end

% now you should define samples of noise and informative signal
noise_samp=input('Enter samples of noise in the form n1:n2\n');
noise=signal(noise_samp,:);
sig=std(noise);

signal_samp=input('Enter samples of signal in the form n1:n2\n');
signal=signal(signal_samp,:);

% define flim
[sp,ffs]=spfft(fs,signal,1);

% flim=input('Enter flim in the form [f1, f2, ..., fn]\n');
% [f0_fft, A0_fft]=choose_fft(sp,ffs,flim);

save(filename, 'noise', 'sig','signal','fs', 'ffs', 'sp','-append')
% figure(2); plot(noise(:,1)); title('noise');
% figure(3); plot(signal(:,1)); title('signal');

%TODO: add sorting signal
