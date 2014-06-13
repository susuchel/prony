function F_sorted = choose_freq(freqs,segSNR,f1,f2)
%CHOOSE_FREQ Chooses frequency estimates within an interval
%   choose_freq(freqs,segSNR,f1,f2) returns an array with two columns
%   First column stands for frequencies, second is for corresponding SNR of
%   the segment which was used to calculate that frequency
%   The array is sorted by SNR in descending order 
%   May be used not only for frequency, but also for damping factors
%   estimate
%
%   Inputs:
%   freqs   - frequency array of size N x M, where N is number of segments,
%           M is number of frequencies found for each segment. Usually N is
%           equal to (N_mod-N_samp+1), M is p/2
%   segSNR  - SNR of each segment, its size is N x 1 that is just a column
%   f1      - lower limit of frequencies
%   f2      - upper limit of frequencies
%   
%   Output:
%   F_sorted    - sorted array with two columns [frequencies, segSNR]
%   
% TODO: change description, because now it deals with observations

%FIXME: if freqs contains all zeros, there will be error, should be fixed
freqs_ch=[]; % new frequency array for chosen frequencies
segSNR_ch=[]; % new array with corresponding SNR

kk=1;

for j=1:size(freqs,3)   % for each observation (third dim) in freqs
    for k=1:size(freqs,2) % for each coloumn in freqs
        for l=1:size(freqs,1) % for each row in freqs
            if ((freqs(l,k,j)<=f2)&&(freqs(l,k,j)>=f1)) % take freq within [f1;f2]
                freqs_ch(kk,1)=freqs(l,k,j);
                segSNR_ch(kk,1)=segSNR(l,j);
                kk=kk+1;
            end
        end
    end
end

F_chosen=[freqs_ch,segSNR_ch]; % chosen data
if ~isempty(F_chosen)
    [~,IX]=sort(F_chosen(:,2),'descend'); % sort by SNR in descending order 
    F_sorted=F_chosen(IX,:); % sorted data
else
    F_sorted=[0 0];
end
end

