function [freqs,alphas] = SProny(signal, fs, p, N_samp,fig_mode, hist_mode, filename)
%SProny Standard Prony procedure 
%   Calculates frequencies and damping factors of a signal
%   Uses least squares procedure for a signal without noise
%   developed 31-05-2014
%   
%   Inputs:
%   signal      - signal to analyze, may be a martix, each column is an 
%               implementation, number of raws = length of the signal 
%   fs          - sampling frequency
%   p           - model order, expected number of exponentials (for real signals 
%               should be 2*n, where n - number of sinusoidals)
%   N_samp      - number of samples to use for estimation
%   fig_mode    - to plot or not to plot roots of polynomial, 
%               if fig_mode~=0, then function plots roots on figure(10*i_impl+1)
%               roots are plotted for each segment, number of roots is equal to p
%               all roots of all segments are plotted on one figure 
%   hist_mode   - to plot or not to plot histogram of frequencies
%               if hist_mode~=0, then function plots freqs of each impl on 
%               figure(10*i_impl+2)
%   filename    - name of file to save freqs, if filename=0, then not saved 
%   
%   Outputs:
%   freqs   - estimated frequencies, size of freqs is L x N x M, 
%   where L=N_samp stands for number of segments, N=p/2 is number of 
%   frquencies, M=size(signal,2) is number of implementations  
%   alphas  - estimates damping factors
%
%   The function estimates frequencies within a segment of N_samp length,
%   then shifts the segment by 1 sample and repeat estimation. Total
%   number of estimates will be (length(signal)-N_samp+1).


N_segm=size(signal,1)-N_samp+1; % number of segments
N_impl=size(signal,2); % number of impementations

freqs=zeros(N_segm,p/2,N_impl); % for each segment p/2 estimates 
                               % should be obtained
alphas=zeros(N_segm,p/2,N_impl); % for each segment p/2 estimates 
                               % should be obtained
                            

                                      
for i_impl=1:N_impl
    
    if fig_mode~=0  %creates figure to plot roots 
        figure(10*i_impl+1);
        hold on
        x=-2:0.00001:2; y=0*x; plot(x,y) % axes
        y=-2:0.00001:2; x=0*y; plot(x,y) 
        % lines, that are borders of a segment where should be frequency 
        % x=-2:0.00001:2; y=2*pi*f1*x/Fs; plot(x,y)
        % x=-2:0.00001:2; y=2*pi*f2*x/Fs; plot(x,y) 
        x=-1:0.00001:1; y=sqrt(1-x.^2); plot(x,y) % unit circle
        y=-sqrt(1-x.^2); plot(x,y)
        xlim([-2 2]); ylim([-2 2]);
    end
    
    for i_segm=1:N_segm    
        S=zeros(N_samp,1);  %segment to analyze 
        for k=1:N_samp      %loop to obtain segment S
            S(k,1)=signal(k+i_segm-1,i_impl); 
        end
        %SNR_S(i_segm,1)=10*log10(sum(S.*S)/(sig*sig*N_samp));
    
    
        % prony procedure starts
        % toeplitz matrix
        Xpf=toeplitz(S(p:N_samp-1),S(p:-1:1)); 

        % toeplitz vector
        xpf=(S(p+1:N_samp));

        % linear prediction coefficients 
        apf=-pinv(Xpf)*xpf;

        % roots of characteristic polynomial 
        A=[1, transpose(apf(1:p))];
        r=roots(A);

        if fig_mode~=0              %plot roots
            plot(r, 'k.')      % plot roots, 'k' stands for black colour
        end

        kk=1;
        for ii=1:length(r)
           if imag(r(ii))>=0     
                if atan(imag(r(ii))/real(r(ii)))/(2*pi*1/fs)>=0
                    freqs(i_segm,kk,i_impl)=atan(imag(r(ii))/real(r(ii)))/(2*pi*1/fs);
                    alphas(i_segm,kk,i_impl)=log(sqrt(imag(r(ii))*imag(r(ii))+real(r(ii))*real(r(ii))))/(1/fs);
                else
                    freqs(i_segm,kk,i_impl)=fs/2+atan(imag(r(ii))/real(r(ii)))/(2*pi*1/fs);
                    alphas(i_segm,kk,i_impl)=log(sqrt(imag(r(ii))*imag(r(ii))+real(r(ii))*real(r(ii))))/(1/fs);
                end
           kk=kk+1;
           end
        end
    end
    
    if hist_mode~=0 
        figure(10*i_impl+2);
        hist(freqs(:,:,i_impl), 1:50:50000);
    end
    
end

if filename~=0
    for i_impl=1:N_impl
    xlswrite(filename, freqs(:,:,i_impl),i_impl); % saves each implem in a sheet
    end
end

end
    

