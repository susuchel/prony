function [mod,sig] = model1 (fs, f0, alfa, N_mod,SNR, N_impl)
%model1 Generates a model signal 
%   model1 (fs, f0, alfa, N_mod,SNR, N_impl) generates a model 
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
%   N_impl - number of model implementations, works if SNR!=0
%   
%   In this model1 magnitudes of the components are same and equal to 0.0249
%   Also initial phases of the components are same and equal to zero
%
%   Each column of the result is an implementation, number of rows=number
%   of samples

freq_op=f0;
a_op=alfa;
n_freq=length(freq_op); % a number of components
t_mod=(0:1/fs:(N_mod-1)/fs); % descrete time

for k=1:n_freq 
    A_mod(1,k)=0.0249; % magnitudes of all components are same 
    fi_mod(1,k)=0; % initial phases of components are same and equal to zero
end

mods=0; % corresponds to the model without noise
for k=1:n_freq
    mods=mods+A_mod(1,k)*(exp(-a_op(k)*t_mod).*cos(2*pi*freq_op(k)*t_mod+fi_mod(1,k)));
end

mod=zeros(N_impl,N_mod); % corresponds to the model with noise
if (SNR==0) 
    mod=mods;
else
    for k=1:N_impl
        sig=sqrt((sum(mods.^2))/(N_mod*10^(SNR/10))); % STD of noise, depending on SNR
        mod(k,:)=mods(1,:)+sig*randn(1,N_mod);
    end
end
 mod=mod'; % 
end