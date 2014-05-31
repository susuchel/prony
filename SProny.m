function [ output_args ] = Untitled3(signal, fs, p, N_samp, )
%SProny Standard Prony procedure 
%   Calculates frequencies and damping factors of a signal
%   Uses least squares procedure for a signal without noise
%   developed 31-05-2014
%   
%   Inputs:
%   signal  - signal to analyze
%   fs      - sampling frequency
%   p       - model order, expected number of exponentials (for real signals 
%           should be 2*n, where n - number of sinusoidals)
%   N_samp  - number of samples to use for estimation
%   
%   Outputs:
%
%   The function estimates frequencies within the interval of N_samp length
%   then shifts the interval by 1 sample and repeat estimation. Total
%   number of estimates will be (length(signal)-N_samp+1).


end

%%%%%%%
Fs=fs; % частота дискретизации сигнала
%%%%%%%%%%%%
qwe1 = mod0; % исследуемый сигнал
%%%%%%%%%%%%
%%%%%%
p=4; % реальное или предполагаемое число экспонент p=2*n, n - число синусоид
%%%%%%%
N=64; % число анализируемых отсчетов 
%%%%%%                        
I0=1; % номер начала анализируемого отрезка
%%%%%%

%-------------------------------------------------------------------------
% Выбираем для анализа отрезок сигнала 
%-------------------------------------------------------------------------

S=zeros(N,1);           % это будет рассматриваемый отрезок сигнала 
for k=0:(N-1) S(k+1)=qwe1(k+I0); end

%-------------------------------------------------------------------------
% Вычисление коэффициентов линейного предсказания
%-------------------------------------------------------------------------

% тёплицева матрицы данных
Xpf=toeplitz(S(p:N-1),S(p:-1:1)); 

% вектор для матричного уравнения
xpf=(S(p+1:N));

% вычисление коэффициентов линейного предсказания 
apf=-pinv(Xpf)*xpf;

%-------------------------------------------------------------------------
% Вычисление корней характеристических полиномов 
%-------------------------------------------------------------------------

% характеристический полином и его корни 
A=[1, transpose(apf(1:p))];
r=roots(A);

figure(3)           % строим корни А

hold on

x=-2:0.00001:2; y=0*x; plot(x,y)
y=-2:0.00001:2; x=0*y; plot(x,y) % наносим на график координатные оси
 
% x=-2:0.00001:2; y=2*pi*20000*x/Fs; plot(x,y)
% x=-2:0.00001:2; y=2*pi*40000*x/Fs; plot(x,y) % прямые, связанные с зоной-сектором,
% %в котором находится искомая частота
 
 x=-1:0.00001:1; y=sqrt(1-x.^2); plot(x,y)
 y=-sqrt(1-x.^2); plot(x,y) % наносим на график единичную окружность
 
plot(r, 'k+')      % наносим на график корни полинома вперед A
xlim([-2 2]); ylim([-2 2]); 
 
% % буква соответствует цвету, символ соответствует символическому обозначению 
% % точки на графике; r - красный, m - фиолетовый, y - желтый, g - зеленый, k
% % - черный, ...
 xlim([-2 2]); ylim([-2 2]); % границы рисунка

%-------------------------------------------------------------------------
% Вычисление частоты и коэффициента затухания
%-------------------------------------------------------------------------

kk=1;
freq=0;
demp=0;
for ii=1:length(r)
       if imag(r(ii))>=0     
        if atan(imag(r(ii))/real(r(ii)))/(2*pi*1/Fs)>=0
        freq(kk)=atan(imag(r(ii))/real(r(ii)))/(2*pi*1/Fs);
        demp(kk)=log(sqrt(imag(r(ii))*imag(r(ii))+real(r(ii))*real(r(ii))))/(1/Fs);
        else
        freq(kk)=Fs/2+atan(imag(r(ii))/real(r(ii)))/(2*pi*1/Fs);
        demp(kk)=log(sqrt(imag(r(ii))*imag(r(ii))+real(r(ii))*real(r(ii))))/(1/Fs);
        end
        kk=kk+1;
       end
end
    
% вывод результата в файл
% %dlmwrite('results.txt', G, '-append', 'delimiter',  '\t', 'precision', '%.6f', 'newline', 'pc')
 xlswrite('freq.xls', freq');
 freq'
 demp'
% %type results.txt
% 
