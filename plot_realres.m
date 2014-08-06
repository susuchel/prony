close all
tic;

SNRlevel=5;%[0,5,10];
svd_mode=1;%[0,1];
N_obs=10;%[1,3];
eps=[0.1, 0.05];%[0.1, 0.05, 0.01];
N_samp=300;
len=625;%1000;
f0=16542;

for z=1:18  % for each factor level
    qwe=RESopt2{z,2};

    for t=2  % for each frequency limit
        T=qwe{t,2};
        str_f0=num2str(f0);
        for k=1:length(N_obs)
            for l=1:length(SNRlevel)
                for i=1:length(svd_mode)
                    Fx=[];
                    Fest=[];
                    Fmean=[];
                    Ffft=[];
                    str_leg={};                
                    for j=1:length(eps)
                        good=(T.SNRlevel==SNRlevel(l)&T.eps==eps(j)&T.svd_mode==svd_mode(i)&T.N_obs==N_obs(k)&T.len==len);
                        Fx(:,j)=T.p(good);
                        Fest(:,j)=T.df_est(good);
                        %Fmean(:,j)=T.df_mean(good);
                        Ffft(:,j)=abs(T.f0_fft(good)-f0)*100/f0;
                        str_leg{j}=sprintf('eps = %g',eps(j));
                    end
                    idfig=z*10000+t*100+i+length(svd_mode)*(l-1)+length(svd_mode)*length(SNRlevel)*(k-1);
                    figure(idfig);
                    strnamefig=sprintf('%d.fig',idfig);
                    strnamebmp=sprintf('%d.bmp',idfig);
                    hold on
                    plot(Fx,Fest);
                    %plot(Fx,Fmean,':');
                    plot(Fx,Ffft,':k');
                    plot(64:144,0.5*ones(length(64:144),1),'k');
                    %legend(str_leg{1},str_leg{2},str_leg{3});
                    legend(str_leg{1},str_leg{2});
                    title(['f0 =',str_f0, ', Nobs = ',num2str(N_obs(k)),', SNRlevel = ',num2str(SNRlevel(l)), ', svdmode= ',num2str(svd_mode(i))]);
                    xlabel('p');
                    ylabel('df, %');
                    ylim([0 2]);
                    saveas(gcf,['C:\Users\Oleg\Documents\GitHub\prony\Figures_thesis\FIG_RES\',strnamefig]);
                saveas(gcf,['C:\Users\Oleg\Documents\GitHub\prony\Figures_thesis\FIG_RES\',strnamebmp]);
                end
            end
        end
    end
end
toc;