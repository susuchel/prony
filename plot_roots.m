function plot_roots( fs,f0, alfa ,flag,varargin)
%plot_roots To plot correct roots in Prony's method
%   TODO: add detailed explanation of inputs

hold on

% plot axes
x=-2:0.00001:2; y=0*x; plot(x,y) 
y=-2:0.00001:2; x=0*y; plot(x,y) 

% plot unit circle
x=-1:0.00001:1; 
y=sqrt(1-x.^2); plot(x,y) 
y=-sqrt(1-x.^2); plot(x,y)
xlim([-2 2]); ylim([-2 2]);

% plot lines passing through the correct roots 
x=-2:0.00001:2;
for k=1:length(f0)
    y=tan(2*pi*f0(k)/fs)*x; plot(x,y)
end

% plot roots of forward ('A') or backward ('B') plynomial in Prony's method
% if flag is 'AB', then plots A and B to subplots
switch flag
    
    case 'A'
        for k=1:length(f0)
            z_A(k)=exp(-alfa(k)/fs+1i*2*pi*f0(k)/fs);
            z_A(length(f0)+k)=exp(-alfa(k)/fs-1i*2*pi*f0(k)/fs);
        end

        if (nargin>4) % user can change LineSpec in plot
            switch nargin
                case 5  
                    plot(z_A,varargin{1})
                case 6
                    plot(z_A,varargin{1},varargin{2})
                case 7
                    plot(z_A,varargin{1},varargin{2},varargin{3})
            end
        else
            plot(z_A,'k^', 'MarkerSize', 12) % default plot settings
        end
        
    case 'B'
        for k=1:length(f0)
            z_B(k)=exp(alfa(k)/fs+1i*2*pi*f0(k)/fs);
            z_B(length(f0)+k)=exp(alfa(k)/fs-1i*2*pi*f0(k)/fs);
        end

        if (nargin>4)
            switch nargin
                case 5  
                    plot(z_B,varargin{1})
                case 6
                    plot(z_B,varargin{1},varargin{2})
                case 7
                    plot(z_B,varargin{1},varargin{2},varargin{3})
            end
        else
            plot(z_B,'k^', 'MarkerSize', 12)
        end
        
    case 'AB'
        for k=1:length(f0)
            z_A(k)=exp(-alfa(k)/fs+1i*2*pi*f0(k)/fs);
            z_A(length(f0)+k)=exp(-alfa(k)/fs-1i*2*pi*f0(k)/fs);
            z_B(k)=exp(alfa(k)/fs+1i*2*pi*f0(k)/fs);
            z_B(length(f0)+k)=exp(alfa(k)/fs-1i*2*pi*f0(k)/fs);
        end
        
        if (nargin>4)
            switch nargin
                case 5  
                    subplot(1,2,1); plot(z_A,varargin{1})
                    subplot(1,2,2); plot(z_B,varargin{1})
                case 6
                    subplot(1,2,1); plot(z_A,varargin{1},varargin{2})
                    subplot(1,2,2); plot(z_B,varargin{1},varargin{2})
                case 7
                    subplot(1,2,1); plot(z_A,varargin{1},varargin{2},varargin{3})
                    subplot(1,2,2); plot(z_B,varargin{1},varargin{2},varargin{3})
            end
        else
            subplot(1,2,1); plot(z_A,'k^', 'MarkerSize', 12)
            subplot(1,2,2); plot(z_B,'k^', 'MarkerSize', 12)
        end
        
    otherwise
    error('In plot_roots(f0,alfa,flag,varargin): flag must be ''A'', ''B'' or ''AB''')
end
end
   
        
         
       