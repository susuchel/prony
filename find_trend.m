close all
tic;

% choose parameters

SNRlevel=5;%[0,5,10];
svd_mode=1;%[0,1];
N_obs=10;%[1,3];
%eps=0.05;%[0.1, 0.05, 0.01];
N_samp=300;
len=625;
%f0=16542;
%p=120;

flim=[13000,15000;
      16000,17000;
      17000,18500;
      26000,27000;
      29000,30000;
      32000,33000;
      36500,38000;
      44000,45000;
      45200,46000;
      57000,59000];

RESest=zeros(10,20);
RESest(:,1:2)=flim;

RESmean=zeros(10,20);
RESmean(:,1:2)=flim;

for p=[120 144]
for eps=[0.1 0.05]

for z=1:18 % for each factor level
    qwe=RESopt2{z,2};

    for t=1:10  % for each frequency limit
        T=qwe{t,2};
        
        good=(T.SNRlevel==SNRlevel&T.eps==eps&T.svd_mode==svd_mode&T.N_obs==N_obs&T.len==len&T.p==p);
        RESest(t,2+z)=T.f_est(good);
        RESmean(t,2+z)=T.f_mean(good);
%                     figure(z*10000+t*100+i+length(svd_mode)*(l-1)+length(svd_mode)*length(SNRlevel)*(k-1));
%                     hold on
%                     plot(Fx,Fest);
%                     %plot(Fx,Fmean,':');
%                     plot(Fx,Ffft,':k');
%                     plot(64:144,0.5*ones(length(64:144),1),'k');
%                     %legend(str_leg{1},str_leg{2},str_leg{3});
%                     legend(str_leg{1},str_leg{2});
%                     title(['f0 =',str_f0, ', Nobs = ',num2str(N_obs(k)),', SNRlevel = ',num2str(SNRlevel(l)), ', svdmode= ',num2str(svd_mode(i))]);
%                     xlabel('p');
%                     ylabel('df, %');
%                     ylim([0 2]);
     
    end
end
toc;

R1=zeros(size(RESest')); R1(1:2,:)=flim'; for i=1:10 R1(3:18,i)=(RESest(i,3:18)'-RESest(i,3))*100/RESest(i,3); end
if eps==0.1
    if p==120
        for i=1:10 figure(i*10+1); plot(R1(3:18,i)); title(['eps = ', num2str(eps), ', p = ', num2str(p), ', ', num2str(R1(1,i)),'  ', num2str(R1(2,i))]); end
    else
        for i=1:10 figure(i*10+2); plot(R1(3:18,i)); title(['eps = ', num2str(eps), ', p = ', num2str(p), ', ', num2str(R1(1,i)),'  ', num2str(R1(2,i))]); end
    end
else     
   if p==120
        for i=1:10 figure(i*10+3); plot(R1(3:18,i)); title(['eps = ', num2str(eps), ', p = ', num2str(p), ', ', num2str(R1(1,i)),'  ', num2str(R1(2,i))]); end
    else
        for i=1:10 figure(i*10+4); plot(R1(3:18,i)); title(['eps = ', num2str(eps), ', p = ', num2str(p), ', ', num2str(R1(1,i)),'  ', num2str(R1(2,i))]); end
    end
end
%R2=zeros(size(RESest')); R2(1:2,:)=flim'; for i=1:10 R2(3:18,i)=(RESmean(i,3:18)'-RESmean(i,3))/RESest(i,3); end
end
end
