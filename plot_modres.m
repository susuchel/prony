close all
tic;
qwe=RES_MOD2;
SNR=[10,5];
N_obs=10; %[1,3,10];
eps=[0.1, 0.05, 0.01];
N_samp=[250,300,400];

alfa=[610, 2000, 1060, 1900, 1480, 870, 3200];
for t=1:7
    T=qwe{t,2};
    str_f0=num2str(table2array(T(1,'f0')));
    for k=1:length(SNR)
        for l=1:length(N_obs)
            for i=1:length(eps)
                Fx=[];
                Fest=[];
                Fmean=[];
                Ffft=[];
                str_leg={};                
                for j=1:length(N_samp)
                    good=(T.SNR==SNR(k)&T.N_samp==N_samp(j)&T.eps==eps(i)&T.N_obs==N_obs(l));
                    Fx(:,j)=T.p(good);
                    Fest(:,j)=T.df_est(good);
                    %Fmean(:,j)=T.df_mean(good);
                    Ffft(:,j)=T.df_fft(good);
                    str_leg{j}=sprintf('Nsamp = %d',N_samp(j));
                end
                idfig=t*100+i+length(eps)*(l-1)+length(eps)*length(N_obs)*(k-1);
                figure(idfig);
                strnamefig=sprintf('%d.fig',idfig);
                strnamebmp=sprintf('%d.bmp',idfig);
                hold on
                plot(Fx,Fest);
                %plot(Fx,Fmean,':');
                plot(Fx,Ffft,':k');
                plot(64:120,0.5*ones(length(64:120),1),'k');
                legend(str_leg{1},str_leg{2},str_leg{3});
                title(['f0 =',str_f0,', alfa = ', num2str(alfa(t)),', SNR = ',num2str(SNR(k)), ', Nobs = ',num2str(N_obs(l)),', eps= ',num2str(eps(i))]);
                xlabel('p');
                ylabel('df, %');
                ylim([0 2]);
                saveas(gcf,['C:\Users\Oleg\Documents\GitHub\prony\Figures_thesis\FIG_MOD\',strnamefig]);
                saveas(gcf,['C:\Users\Oleg\Documents\GitHub\prony\Figures_thesis\FIG_MOD\',strnamebmp]);
            end
        end
    end
end
toc;