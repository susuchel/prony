function [mod,sig] = model1 (fs, f0, alfa, N_mod,SNR, N_obs)
%model1 Generates a model signal 
%   model1 (fs, f0, alfa, N_mod,SNR, N_obs) generates a model 
%   that is a sum of sinusoidal components. Output is N_mod samples.
%   
%   Input arguments:
%   fs      - sampling frequency
%   f0      - a vector of frequecies of each sinusoidal component
%   alfa    - a vector of damping factors of each component, should be
%   positive numbers in case there is damping
%   N_mod   - number of samples
%   SNR     - signal-to-noise ratio that is defined as
%   SNR=10*log10(sum(s[n]^2)/(N*sig^2)), s[n] - samples of a signal without
%   noise, sig - std of noise
%   if SNR = 0 then model without noise will be generated
%   N_impl - number of obseravtions, works if SNR!=0
%
%   If some of the inputs are [], then default values will be used
%   
%   In this model1 magnitudes of the components are same and equal to 0.0249
%   Also initial phases of the components are same and equal to zero
%
%   Each column of the result is an observation, number of rows=number
%   of samples
%   TODO: change description because now sig is a vector

% define default values, works if input argument is empty
if isempty(fs)
    fs=500000;
end

if isempty(f0)
    f0=[13853, 17690, 27690, 29600, 32600, 35421, 44310];
    %   Amplitudes may be used later:
    %   A = [0.0249, 0.0742, 0.03, 0.0347, 0.0522, 0.0216, 0.2487];
end

if isempty(alfa)
    alfa=[610, 2000, 1060, 1900, 1480, 870, 3200];
end

if isempty(N_mod)
    N_mod=2500;
end

freq_op=f0;
a_op=alfa;
n_freq=length(freq_op); % a number of components
t_mod=(0:1/fs:(N_mod-1)/fs); % descrete time

A_mod=zeros(1,n_freq);  
fi_mod=zeros(1,n_freq);
for k=1:n_freq 
    A_mod(1,k)=0.00249; % magnitudes of all components are same 
    fi_mod(1,k)=0; % initial phases of components are same and equal to zero
    % further idea is to use random A_mod and fi_mod
end

mods=0; % corresponds to the model without noise
for k=1:n_freq
    mods=mods+A_mod(1,k)*(exp(-a_op(k)*t_mod).*cos(2*pi*freq_op(k)*t_mod+fi_mod(1,k)));
end

mod=zeros(N_obs,N_mod); % corresponds to the model with noise
sig=zeros(N_obs,1);
sig1=0;
if (SNR==0) 
    mod=mods;
else
    sig1=sqrt((sum(mods.^2))/(N_mod*10^(SNR/10))); % STD of noise, depending on SNR
    for k=1:N_obs
        mod(k,:)=mods(1,:)+sig1*randn(1,N_mod);
        sig(k,1)=sig1;
    end
end
 mod=mod'; % 
end