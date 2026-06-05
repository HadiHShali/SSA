function plot_n(varargin)

    t=varargin{1};
    legendlabel=varargin{end};
    numplots=length(legendlabel); %the legend has the total number of plots (data and residuals - with residuals calculated here)
    n_vec_in=nargin-2; %number of input vectors to plot
    %(nargin has first argument for time vector and last for legend, the
    %rest are input vectors to plot

    if numplots==2
        symcolr={'c' 'r--'};
    else
        symcolr={'c' 'r' 'b' 'k' 'm'};
    end

    hold on
    for ydata_n=1:n_vec_in
        y(:,ydata_n)=varargin{ydata_n+1};
        plot(t,y(:,ydata_n),symcolr{ydata_n})
    end    
    xlabel('time (days)')
    ylabel('position')

    numres=numplots-n_vec_in;
    scl=10;
    offset=-3;
    for resno=1:numres
        res=y(:,1)-y(:,resno+1);
        maxres=max(abs(res));
        plot(t,scl*(y(:,1)-y(:,resno+1))+offset,symcolr{resno+2})
    end

    legend(legendlabel,'Location','southeast')

end