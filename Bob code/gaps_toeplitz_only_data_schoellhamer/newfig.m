function [figno fig_handle]=newfig(figno,figname)
    figno=figno+1;
    fig_handle=figure('Name',figname);
end