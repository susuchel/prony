function plot_results(T,C,var1, var2)
%plot_results(T,C,var1, var2) Summary of this function goes here
%   Detailed explanation goes here

Var1=table2array(T(:,var1));
Var2=table2array(T(:,var2));
plot(Var1,Var2, 'go');

end

