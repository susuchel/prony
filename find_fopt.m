function fopt = find_fopt( f0 )
%FIND_OPT Finds optimal sampling rate in Prony's Method 
%   find_fopt( f0 ) For given number of frequency f0 it calculates
%   corresonding sampling rate fopt
%   For more details see: 
%   Bushuev, Ibryaeva (2012, July). DOI 10.1109/TSP.2012.6256374 
%

f0=sort(f0);    % f0 must be sorted in acsending order

Fnum=0;
Fdenum=0;

for k=1:length(f0)
    Fnum=Fnum+f0(k)*f0(k);          % numerator
    Fdenum=Fdenum+(2*k-1)*f0(k);    % denumerator
end

fopt=4*length(f0)*Fnum/Fdenum;

if (fopt<=2*f0(end))    % to be useful, fopt should be greater than 2*max(f0)
    warning('fopt/2 is less than input frequency %g Hz', f0(end));
end
end

