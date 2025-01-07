clear
close all

%format longG %decimal to 18 sig digits, then sci notation
format shortG %decimal to 7 sig digits, then sci notation

figno=0;

%process hector file from Hadi
filename='1LSU_1_SSA.dat'; %input data has gaps (time not continuous)
[data datasize]=read_matrix(filename);
%big natural gap in data from 4847 (last day) to 5085 (first day back, 238-1 missing days).
n_data=size(data,1); %this is number data points (not same as number days from start to finish due to gaps)

strtdata=1;%1000;
enddata=n_data;%4600;

%x_in=data(:,3); %ordered but not continuous
x_in=data(strtdata:enddata,3); %ordered but not continuous
%tyfrac=data(:,2); %get time in fraction of year, this good for plotting
tyfrac=data(strtdata:enddata,2); %get time in fraction of year, this good for plotting
%get time series, consists of both sample time of data and data (not like seismograms that are continuous with no gaps
%tmjd=data(:,1); %get time in modified julian days, this good for indexing
tmjd=data(strtdata:enddata,1); %get time in modified julian days, this good for indexing
dataindx=tmjd-tmjd(1)+1; %this has indices of days with data (will line up with continuous t later)

%remove mean and normalize here (as opposed to after filling in gaps with NaNs)
[x_std x_mean]=std(x_in);
x=(x_in-x_mean)/x_std; %this is data for ssa

t=[tmjd(1):tmjd(end)]-tmjd(1)+1; %make a continuous time vector that starts counting from 1, indices in dataindx are the elemenets with data
n_cont_time=length(t); %length (number days) of continuous time vector that is ≥n_data (= if no gaps, else n_cont_time>n_data)
x_with_NaN=nan(n_cont_time,1); %make vector of NaNs of continuous time length
x_with_NaN(dataindx)=x; %put data into continuous series where observations have values, missing data remain NaNs
%x_with_NaN(4800+[0:1273])=NaN; %test - make ts with big gap longer than embedding dimension (0-1273=1274 points, is good, at 1274/1275 points get missing elements in PCs) - get missing PC terms

x = x_with_NaN;        %this has mean removed and normalized to standard deviation, did not have to process the NaNs (eg. mean omitting NaNs) to get this by doing before inserting NaNs.

%-----------------------------------------------------------

figname='gappy input time series, mean removed and normalized, NaNs for missing data';
figno=newfig(figno,figname);
plot(t,x,'bo-');
xlabel('time (days)')
ylabel('position')
closefig(figname)

%disp('find gaps (NaNs)')
[NaNindx logical]=find(isnan(x));
num_nans=length(NaNindx); %get number nans
frac_nans=num_nans/n_cont_time; %fraction of NaNs for one step of how bad NaNs are (problem not lots individual NaNs, but long sequences)

n_data_ck=n_cont_time-num_nans;
diff_ck(n_data,n_data_ck,"check num data vs num gappy data in continuous");

%m=fix(n/2); %m is number columns in traj matrix, not number points in each column
%m=fix(n_data/2); %can only have as many RC as [fix(number data/2)]. m≤half number data (and continuous length≥ndata)
%this is (generally) max embedding dimension (max number RCs and shortest
%sections of time series)
m=1278;%3.5 years
%m=1278; %set embedding dimension to 3.5 years (rounded to integer), gives number of columns in trajectory matrix and number RCs we will obtain
l=n_cont_time-m+1; %number rows in trajectory matrix (for m=n/2, with n even, is one taller than wide, with n odd is 2 taller than wide), as m is smaller, l is longer

colind=1:m; %colmn index vector
rowind=0:l-1; %row index vector
trajmatind=colind+rowind'; %vectorized build of trajectory matris
trajmat=x(trajmatind);

OutStrucNaNs=SSABasicCleanUCLA(trajmat); %locally defined below - fixed for missing data

%how to pick this automatically
% how to pick this automatically????\
maxeig=29;
%maxeig=100;
%maxeig=2;
go=1;

while go %0 is false, anything else true

    figname=strcat("eigenvalues, max eig=",num2str(maxeig));
    figno=newfig(figno,figname);
    semilogy(OutStrucNaNs.LAMBDA,'b+-')
    hold on
    semilogy(OutStrucNaNs.LAMBDA(1:maxeig),'ro-')
    closefig(figname)
    
    RCunscaled=OutStrucNaNs.RCunscaled;
    RCscaled=OutStrucNaNs.RCscaled;
    RCunscaledNaNs=InsertNaNs(RCunscaled,NaNindx); %puts NaNs into the RCs when there are NaNs in the input data
    RCscaledNaNs=InsertNaNs(RCscaled,NaNindx); %puts NaNs into the RCs when there are NaNs in the input data
    recons_unscaledNaNs=sum(RCunscaledNaNs,2);
    recons_scaledNaNs=sum(RCscaledNaNs,2);
    diff_ck(x,recons_unscaledNaNs,"x-recons_unscaledNaNs")
    diff_ck(x,recons_scaledNaNs,"x-recons_scaledNaNs")
    %x_reconsNaNs=reconsNaNs*x_std_omitnan+x_mean_omitnan;
    %diff_ck(x_in,x_reconsNaNs,"x_in-x_reconsNaNs")
    
    figname=strcat("recon_unscaled gappy: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
        ", M=", num2str(m), ", L=", num2str(l), ", max eig=", num2str(maxeig));
    figno=newfig(figno,figname);
    plot(t,x,'b')
    hold on
    plot(t,recons_unscaledNaNs,'r--')
    rc_unscaled=sum(RCunscaledNaNs(:,[1:maxeig]),2);
    rc_unscaled(NaNindx)=NaN;
    plot(t,rc_unscaled,'k')
    xlabel('time (days)')
    ylabel('position')
    legend({'in' 'full RC' 'gappy sig model'},'Location','southeast')
    closefig(figname)
    
    figname=strcat("recon_scaled gappy: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
        ", M=", num2str(m), ", L=", num2str(l), ", max eig=", num2str(maxeig));
    figno=newfig(figno,figname);
    plot(t,x,'b')
    hold on
    plot(t,recons_scaledNaNs,'r--')
    rc_scaled=sum(RCscaledNaNs(:,[1:maxeig]),2);
    rc_scaled(NaNindx)=NaN;
    plot(t,rc_scaled,'k')
    xlabel('time (days)')
    ylabel('position')
    legend({'in' 'full RC' 'gappy sig model'},'Location','southeast')
    closefig(figname)

    figure
    figname="reconstructions";
    figno=newfig(figno,figname);
    rc_scaled=sum(RCscaledNaNs(:,[1:maxeig]),2);
    rc_scaled(NaNindx)=NaN;
    plot(t,rc_scaled,'b')
    hold on
    rc_unscaled=sum(RCunscaledNaNs(:,[1:maxeig]),2);
    rc_unscaled(NaNindx)=NaN;
    plot(t,rc_unscaled,'r--')
    xlabel('time (days)')
    ylabel('position')
    legend({'scaled' 'unscaled'},'Location','southeast')
    closefig(figname)
    
    diff_ck(rc_scaled,rc_unscaled,"rc_scaled-rc_unscaled")

    figname=strcat("recon_unscaled gappy TS residual NOISE, max eig+1=", num2str(maxeig+1)," to end");
    figno=newfig(figno,figname);
    ruidoest=sum(RCunscaledNaNs(:,[maxeig+1:end]),2);
    plot(t,ruidoest,'r')
    legend({figname},'Location','southeast')
    closefig(figname)

    figname=strcat("recon_scaled gappy TS residual NOISE, max eig+1=", num2str(maxeig+1)," to end");
    figno=newfig(figno,figname);
    ruidoest=sum(RCscaledNaNs(:,[maxeig+1:end]),2);
    plot(t,ruidoest,'r')
    legend({figname},'Location','southeast')
    closefig(figname)



    %go=input(strcat("current maxeig=", num2str(maxeig)," enter new maxeig, <CR> for quit "));
    go=0;
end

