function [freqs,alphas,segSNR] = MProny(signal, fs, sig, p, q, N_samp,eps,fig_mode, hist_mode, filename)
%MProny Modified Prony procedure 
%   Calculates frequencies and damping factors of a signal
%   Uses least squares procedure for a signal
%   Uses several tricks to improve performance
%   developed 04-06-2014
%   
%   Inputs:
%   signal      - signal to analyze, may be a martix, each column is an 
%               observation, number of raws = length of the signal 
%   fs          - sampling frequency
%   sig         - std of noise of a signal, may be different for each
%               observation, but now it is scalar (to be changed when 
%               analyzing real signals)
%   p           - model order, higher than expected number of exponentials
%   q           - expected number of exponentials (for real signals 
%               should be 2*n, where n - number of sinusoidals)
%   N_samp      - number of samples to use for estimation
%   eps         - epsilon, difference btw roots of 'A' and 'B'
%   fig_mode    - to plot or not to plot roots of polynomial, 
%               if fig_mode~=0, then function plots roots on figure(10*i_obs+1)
%               roots are plotted for each segment, number of roots is equal to p
%               all roots of all segments are plotted on one figure 
%   hist_mode   - to plot or not to plot histogram of frequencies
%               if hist_mode~=0, then function plots freqs of each obs on 
%               figure(10*i_obs+2)
%   filename    - name of file to save freqs, if filename=0, then not saved 
%   
%   Outputs:
%   freqs   - estimated frequencies, size of freqs is L x N x M, 
%   where L=N_samp stands for number of segments, N=q/2 is number of 
%   frequencies (may be of different size), M=size(signal,2) is number of 
%   observations  
%   alphas  - estimated damping factors, absolute values (positive)
%   segSNR  - SNR of each segment
%
%   The function estimates frequencies within a segment of N_samp length,
%   then shifts the segment by 1 sample and repeat estimation. Total
%   number of estimates will be (length(signal)-N_samp+1).


N_segm=size(signal,1)-N_samp+1; % number of segments
N_obs=size(signal,2); % number of observations

freqs=zeros(N_segm,p/2,N_obs); % for each segment p/2 estimates 
                               % should be obtained
alphas=zeros(N_segm,p/2,N_obs); % for each segment p/2 estimates 
                               % should be obtained
segSNR=zeros(N_segm,N_obs);  % SNR of each segment                          

                                      
for i_obs=1:N_obs
    
    if fig_mode~=0  %creates figure to plot roots 
        figure(10*i_obs+1);
        
        % to plot roots of A
        subplot(1,2,1); 
        hold on
        % plot axis
        x=-2:0.00001:2; y=0*x; plot(x,y) 
        y=-2:0.00001:2; x=0*y; plot(x,y) 
        % plot unit circle
        x=-1:0.00001:1; 
        y=sqrt(1-x.^2); plot(x,y) 
        y=-sqrt(1-x.^2); plot(x,y)
        xlim([-2 2]); ylim([-2 2]);
                    
        % to plot roots of B
        subplot(1,2,2); 
        hold on
        % plot axis
        x=-2:0.00001:2; y=0*x; plot(x,y) 
        y=-2:0.00001:2; x=0*y; plot(x,y) 
        % plot unit circle
        x=-1:0.00001:1; 
        y=sqrt(1-x.^2); plot(x,y) 
        y=-sqrt(1-x.^2); plot(x,y)
        xlim([-2 2]); ylim([-2 2]);
    end
    
    for i_segm=1:N_segm    
        S=zeros(N_samp,1);  %segment to analyze 
        for k=1:N_samp      %loop to obtain segment S
            S(k,1)=signal(k+i_segm-1,i_obs); 
        end
        segSNR(i_segm,i_obs) = calc_SNR2 (S,sig); % suppose sig is a scalar
        %SNR_S(i_segm,1)=10*log10(sum(S.*S)/(sig*sig*N_samp));
    
    
        % prony procedure starts
        % toeplitz matrix
        Xpf=toeplitz(S(p:N_samp-1),S(p:-1:1)); % forward
        Xpb=toeplitz(S(p+1:N_samp),S(p+1:-1:2)); % backward

        % toeplitz vector
        xpf=(S(p+1:N_samp)); % forward
        xpb=(S(1:N_samp-p)); % backward

        % linear prediction coefficients 
        apf=-pinv(Xpf)*xpf; % forward
        apb=-pinv(Xpb)*xpb; % backward

        % roots of characteristic polynomial 
        A=[1, transpose(apf(1:p))]; % forward
        r=roots(A);
        
        B=[1, ctranspose(apb(p:-1:1))]; % backward
        conjr=roots(B);

        if fig_mode~=0              %plot roots
            subplot(1,2,1);
            plot(real(r), imag(r), 'k.', 'MarkerSize', 5)      
            
            subplot(1,2,2);
            plot(real(conjr), imag(conjr), 'k.', 'MarkerSize', 5) 
        end
        
        
        chosen=0;       % for chosen roots of B (imag>0 and abs>1) 
        implus=0;
        ii=1; iii=1;        

        for k=1:length(conjr)
            if imag(conjr(k))>0
                implus(ii)=conjr(k);  
                if abs(implus(ii))>1      
                    chosen(iii)=implus(ii); 
                    iii=iii+1;
                end
                ii=ii+1;
            end
        end
        
        for kk=1:length(chosen)
            for ii=1:length(r)
                if abs(chosen(kk)*r(ii)-1)<eps   % choose true roots
                    if atan(imag(r(ii)')/real(r(ii)'))/(2*pi*1/fs)>=0
                        freqs(i_segm,kk,i_obs)=atan(imag(r(ii)')/real(r(ii)'))/(2*pi*1/fs);
                        alphas(i_segm,kk,i_obs)=-log(sqrt(imag(r(ii)')*imag(r(ii)')+real(r(ii)')*real(r(ii)')))/(1/fs);
                    else
                        freqs(i_segm,kk,i_obs)=fs/2+atan(imag(r(ii)')/real(r(ii')))/(2*pi*1/fs);
                        alphas(i_segm,kk,i_obs)=-log(sqrt(imag(r(ii)')*imag(r(ii)')+real(r(ii)')*real(r(ii)')))/(1/fs);
                    end
           	
                end
            end
        end
    end
    
    if hist_mode~=0
        figure(10*i_obs+2);
        hist(freqs(:,:,i_obs), 1:50:50000);
    end
end

if filename~=0
    for i_obs=1:N_obs
    xlswrite(filename, freqs(:,:,i_obs),i_obs); % saves each obs in a sheet
    end
end

end
    

