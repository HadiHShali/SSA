clear
close all

%format longG %decimal to 18 sig digits, then sci notation
format shortG %decimal to 7 sig digits, then sci notation

figno=0;

%process hector file from Hadi

filename='1LSU_1_SSA.dat'; %(all real) input data has gaps (mjd, fracyr, posn; time not continuous)
[full_data]=read_matrix(filename);
%big natural gap in data from 4847 (last day) to 5085 (first day back, 238-1 missing days).

n_data=size(full_data,1); %this is number data points (not same as number days from start to finish due to gaps)

%strtdata=1;
%enddata=length(full_data);
data=full_data;

x_in=data(:,3); %ordered but not uniformly spaced, column vector

tyfrac=data(:,2); %time in fraction of year (from input file), this good for plotting, not uniformly spaced, column vector
%get time series, consists of both sample time of data and data (not like seismograms that are continuous with no gaps

tmjd=data(:,1); %get time of data in modified julian days, this good for indexing, not uniformly spaced, column vector
dataindx=tmjd-tmjd(1)+1; %this has indices of days with data (will line up with continuous t later)

%remove mean and normalize here (as opposed to after filling in gaps with NaNs)
[x_std x_mean]=std(x_in);
x_in_mean_rem_norm_w_std=(x_in-x_mean)/x_std; %this is data for ssa, mean removed and normalized by sd

t=([tmjd(1):tmjd(end)]-tmjd(1)+1)'; %make a continuous time vector that starts counting from 1, indices in dataindx are the elemenets with data
%no NaNs in t
n_cont_time=length(t); %length (number days) of continuous time vector that is ≥n_data (= if no gaps, else n_cont_time>n_data)
if n_cont_time==n_data
    disp('No NaNs in input data')
else
    disp(strcat("there are ", num2str(n_cont_time-n_data)," NaNs"))
end

N=n_cont_time;

x_in_with_NaNs=nan(n_cont_time,1); %make copy input vector of continuous time length with NaNs
x_in_with_NaNs(dataindx)=x_in; 

x_in_mean_rem_norm_w_std_with_NaNs=nan(n_cont_time,1);
x_in_mean_rem_norm_w_std_with_NaNs(dataindx)=x_in_mean_rem_norm_w_std; %put data into continuous series where observations have values, missing data remain NaNs

x = x_in_mean_rem_norm_w_std_with_NaNs;%this has mean removed and normalized to standard deviation, did not have to process the NaNs (eg. mean omitting NaNs) to get this by doing before inserting NaNs.

x_in_with_Zeros_mean_rem_norm_w_std=x_in_mean_rem_norm_w_std_with_NaNs;
x_in_with_Zeros_mean_rem_norm_w_std(isnan(x_in_mean_rem_norm_w_std_with_NaNs))=0;
X=x_in_with_Zeros_mean_rem_norm_w_std; %NaNs replaced with zeros (no zeros in real data) so don't have to play with NaNs

%-----------------------------------------------------------

figname='gappy input time series, mean removed and normalized, NaNs for missing data';
[figno fig_handle]=newfig(figno,figname);
plot(t,x,'bo-');
xlabel('time (days)')
ylabel('position')
showfig(figname)

%disp('find gaps (NaNs)')
[NaNindx logical]=find(isnan(x));
num_nans=length(NaNindx); %get number nans
frac_nans=num_nans/n_cont_time; %fraction of NaNs for one step of how bad NaNs are (problem not lots individual NaNs, but long sequences)

n_data_ck=n_cont_time-num_nans;
diff_ck(n_data,n_data_ck,"check num data vs num gappy data in continuous");

%embedding dimension - gives number COLUMNS in traj matrix, not number points in each column
%the embedding dimension (m) sets the "frequency response" (gut feeling - periods much longer go into trends,
%much shorter into oscillations/cycles)

%M=fix(n_data/2); %can only have as many RC as [fix(number data/2)] so m≤half number of data (and continuous length≥ndata)
%line above is (generally) max embedding dimension (max number RCs and shortest sections of time series)

M=1278; %set embedding dimension to 3.5 years (rounded to integer),
%gives number of columns in trajectory matrix and number RCs we will obtain
%M=m;
l=n_cont_time-M+1; %number rows in trajectory matrix (for m=n/2, with n even, is one taller than wide, with n odd is 2 taller than wide), as m is smaller, l is longer

colind=1:M; %colmn index vector
rowind=0:l-1; %row index vector
trajmatind=colind+rowind'; %vectorized build of trajectory matrix
trajmat=NaN(l,M);
trajmat=x(trajmatind);

%do with toeplitz matrix
%EACH C(Jindx) IS INDEPENDENT OF THE OTHERS SO THE LOOP ON J CAN BE PARALLIZED
%this is faster than vectorizing (at least as far as I've tried to vectorizie it)
for J=0:M-1
    Jindx=J+1;
    Cloop(Jindx)=0;
    Nterms=N-J;
    cntZeros=0;
    for I=1:Nterms
        newProdTerm=X(I)*X(I+J);
        if newProdTerm~=0
            Cloop(Jindx)=Cloop(Jindx)+newProdTerm;
        else
            cntZeros=cntZeros+1;
        end
    end
    ActualNterms=Nterms-cntZeros;
    Cloop(Jindx)=Cloop(Jindx)/ActualNterms;
end
figname=strcat("toeplitz diagonal values double loop - mine 1, gappy: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
    ", M=", num2str(M), ", L=", num2str(l));
figno=newfig(figno,figname);
plot(Cloop)

Ctoep=toeplitz(Cloop);

figname=strcat("toeplitz lagged correlation, gappy: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
    ", M=", num2str(M), ", L=", num2str(l));
figno=newfig(figno,figname);
imagesc(Ctoep)

%OutStrucNaNs=SSABasicCleanUCLAtoeplitz(trajmat,Ctoep);
OutStrucNaNs=SSABasicCleanUCLAtoeplitz(trajmat);

%how to pick this automatically
% how to pick this automatically????
maxeig=14;
maxeig=15;
%maxeig=18;
maxeig=29;
%maxeig=50;
%maxeig=100;
%maxeig=m/2;
%maxeig=2;

figname=strcat("eigenvalues, max eig=",num2str(maxeig));
figno=newfig(figno,figname);
semilogy(OutStrucNaNs.LAMBDA,'b+-')
hold on
semilogy(OutStrucNaNs.LAMBDA(1:maxeig),'ro-')
showfig(figname)

RC=OutStrucNaNs.RC; %to do this "properly" (without cut and paste duplication of the 2 processings requires 3d arrays
RCNaNs=InsertNaNs(RC,NaNindx); %puts NaNs into the RCs when there are NaNs in the input data
%above replaces the zeros that are place holders for NaNs with NaNs again
full_recons_NaNs=sum(RCNaNs,2); %reconstruct full time series for 
x_in_recons_NaNs=full_recons_NaNs*x_std+x_mean;%mean and std computed before inserting NaNs,
%unnormalize and reinsert mean to compare with input (also need to process this to get velocity)

figname=strcat("gappy full reconstructions and data");
figno=newfig(figno,figname);
plot_n(t,x,full_recons_NaNs,{'input data','full reconstruction'})
showfig(figname)

figname=strcat("residual");
figno=newfig(figno,figname);
plot_n(t,x-full_recons_NaNs,{'residual'})
showfig(figname)

rc=sum(RCNaNs(:,[1:maxeig]),2);
rc(NaNindx)=NaN;

figname=strcat("part  recon: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
    ", M=", num2str(M), ", L=", num2str(l), ", num eig=", num2str(maxeig));
figno=newfig(figno,figname);
%plot_n(t,x,rc,{'data' 'model '})
plot_n(t,x,rc,{'data' 'model '})
showfig(figname)

figname=strcat("multiple part  recon: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
    ", M=", num2str(M), ", L=", num2str(l), ", num eig=", num2str(maxeig));
figno=newfig(figno,figname);
plot_n(t,x,sum(RCNaNs(:,[1:2]),2),sum(RCNaNs(:,[3:4]),2), ...
    sum(RCNaNs(:,[5:6]),2),sum(RCNaNs(:,[7:maxeig]),2),{'data' '1-2' '3-4' '5-6' strcat('7-', num2str(maxeig))})
showfig(figname)

ruidoestRCNaNs=sum(RCNaNs(:,[maxeig+1:end]),2);

figname=strcat("recon_ gappy TS residual NOISE, max eig+1=", num2str(maxeig+1)," to end");
figno=newfig(figno,figname);
plot(t,ruidoestRCNaNs,'r')
legend({figname},'Location','southeast')
showfig(figname)

%check this out - since I put zeros in the traj matrix in the gaps, one gets zeros back
%in the full reconstruction (that I removed where ever there was a NaN above) .... eventually
figname=strcat("full  cont recon and data");
figno=newfig(figno,figname);
for ntrms=[2:20:M M]
    clf
    plot_n(t,x,sum(RC(:,1:ntrms),2),{'input data',strcat('continuous full reconstruction: ', num2str(ntrms))})
    set(gcf, 'Position', get(0, 'Screensize'));
    title(strcat(figname,", number terms ",num2str(ntrms)))
    drawnow
end
%showfig(figname)

%BUT the low order RC are smoothly continuous through the gaps!! This is nice for fourier analysis of the
%low order RC (that are close to sinusoidal)
%probably want to play with max rc so interpolation is smoother than that in the 7:maxeig curve
figname=strcat("mult part cont recon: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
    ", M=", num2str(M), ", L=", num2str(l));
figno=newfig(figno,figname);
plot_n(t,x,sum(RC(:,[1:2]),2),sum(RC(:,[3:4]),2),sum(RC(:,[5:6]),2),sum(RC(:,[7:maxeig]),2),{'data' '1-2' '3-4' '5-6' '7-end'})
showfig(figname)
