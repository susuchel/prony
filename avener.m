function y = avener( x )
%avener Average energy of a signal   
%   avener(x) finds average energy of a signal x
%   if x is a matrix, avener calculates average energy of each column

y=zeros(1,size(x,2));
for k=1:size(x,2)
    y(1,k)=sum(x(:,k).*x(:,k))/size(x,1);
end


end

