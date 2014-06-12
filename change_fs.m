function [signal_new, fs_new] = change_fs( signal, fs, fopt )
% CHANGE_FS Changes sampling rate of a signal
%   change_fs(signal, fs, fopt) changes the sampling rate of a signal from
%   initial fs to fs_new, that is close to fopt and is in the form of
%   fs/koeff, where koeff - an integer
%   Each initial observation after resampling will result in (koeff) new 
%   observations, so the total number of observations will be N_obs*koeff
%   New observations are put in signal_new in the following way:
%   each new observation is a column and includes N_new samples,
%   first N_obs columns are resampled initial observations started from 
%   first sample, next N_obs columns are resampled initial observations 
%   started from second sample and so on. In other words first initial
%   observation is split into (koeff) new observations, first new
%   obseravation is in first column, second is in (N_obs+1) column, etc.

%   Determine koeff to resample signal
%   Suppose fs/fopt=(A+r), A is an integer and r is remainder over fopt
%   if r > A/(2A+1), then fs/(A+1) is closer to fopt than fs/A
%   if r < A/(2A+1), then one should take fs/A

A=fix(fs/fopt); % integer part of fs/fopt
r=fs/fopt-A;    % the remainder over fopt

if (r>=A/(2*A+1))
    koeff=A+1;
else
    koeff=A;
end
    
fs_new=fs/koeff;

% Do resampling

[N_old, N_obs]=size(signal); % N_obs stands for number of obseravations of 
                            % the initial signal (number of columns), 
                            % N_old - number of samples in each obseravation
                            % (number of rows)
N_new=floor(N_old/koeff);   % N_new is number of samples after resampling
                            % Each observation after resampling will result
                            % in (koeff) new observations, so the total
                            % number of observations will be N_obs*koeff

signal_new=zeros(N_new,koeff*N_obs);

%for k=1:(floor((length(signal)-1)/koeff)+1) 
% New observations are put in signal_new in the following way:
% each new observation is a column and includes N_new samples
% first N_obs columns are resampled observations started from first sample
% following N_obs columns are resampled observations started from second sample
% and so on
for i_obs=1:N_obs
    for j=1:koeff
        for k=1:N_new
            signal_new(k, i_obs+(j-1)*N_obs)=signal(j+koeff*(k-1),i_obs);
        end
    end
end
% TODO: consider 3-dimensional array, or output koeff, otherwise it is not
% clear how many observations from first sample are in output signal_new 
% i.e. N_new x N_obs x koeff
end

