clear
close all

ts=getdata;
m=1278; %set embedding dimension to 3.5 years (rounded to integer),
%gives number of columns in trajectory matrix and number RCs we will obtain
f=0;

[xmean,xstd,xk,ro,lam,xk_]=vssanan(ts,m,f);
xf=rc(xk,ro);
xf_=rc(xk_,ro);

nan_indx=find(isnan(ts))';
RC=(sum(xf)*xstd)+xmean;
RC_=(sum(xf_)*xstd)+xmean;
RC(nan_indx)=nan;
RC_(nan_indx)=nan;

max(abs(ts-RC'))
max(abs(ts-RC_'))

plot(ts,'b')
hold on
plot(RC,'g')
plot(RC_,'r--')

%--------------------------------------------------------------------------

function [xmean,xstd,xk,ro,lam,xk_]=vssanan(ts,m,f)
    %function [xmean,xstd,xk,ro,lam]=vssanan(ts,m,f)
    %
    % purpose: given time series ts and window size m (m*dt), compute
    %          principal components xk, eigenvectors ro, and eigenvalues lam.
    %          see Vautard and Ghil 1989 for notation
    %          this version tries to ignore nans
    %
    % inputs: ts=time series with constant dt
    %         m=window size, in time steps
    %         f=fraction (0<f<=1) of good data points for determining pc's
    %        
    % outputs: xmean=mean of ts
    %          xstd=standard deviation of ts
    %          xk=principal components zero mean time series, row k = pc k
    %          ro=eigenvectors ro(j,k), k=corresponding pc
    %          lam=eigenvalue vector, sorted from greatest to smallest
    %
    % by David Schoellhamer, 5/27/93, 12/99, 1/2000
    %
    % normalize time series x
    %
    igood=find(~isnan(ts));
    xmean=mean(ts(igood));
    xstd=std(ts(igood));
    x=(ts-xmean)/xstd;
    n=length(x);
    %
    % compute autocorrelation function eqn b.2
    %
    % only need the first m values
%    disp('Compute autocorrelation function')
    for j1=1:m
        %summ=0;
        j=j1-1;
        prod=x(1:n-j).*x(j+1:n);
        igood=find(~isnan(prod));
        c(j1)=sum(prod(igood))/(length(igood)-1);
    end
    %
    %form the covariance matrix a, divide by m (eqn b.3)
    %
    %disp('Form covariance matrix a')
    a=toeplitz(c)/m;
    %
    %determine eigenvalues and eigenvectors of a (solve b.3)
    %
 %   disp('Determine eigenvalues and eigenvectors of a');
    [z,eval]=eig(a);
    lam=diag(eval);
    %
    % sort by descending eigenvalues
    %
    %disp('Sort by descending eigenvalues');
    [lam,ilam]=sort(lam);
    lam=flipud(lam);
    ilam=flipud(ilam);
    summ=0;
    kmax=0;
    semilogy(lam(1:min([m 30])),'o')
    nkmax=0;

    nkmax=m-1; %rs
    %nkmax=13; %rs - to make it go fast while testing other stuff (see eigenvalue plot)

    while nkmax==0
        nkmax=input('Number of pcs (0 for keyboard control) = ');
        if nkmax==0
            disp('Enter return command to resume')
            keyboard
        end
    end 
    while kmax<nkmax
        kmax=kmax+1;
        ro(:,kmax)=z(:,ilam(kmax));
        summ=summ+lam(kmax);
        %disp(sprintf('Eigenvalue of pc %d is %g, sum of variance is %g',kmax,lam(kmax),summ))
    end
    %
    % determine principal component time series (eqn b.6)
    %
    xk=nan*ones(kmax,n-m+1);
    xk_=nan*ones(kmax,n-m+1); %rs, principal components not scaled
    
    for k=1:kmax
        %disp(sprintf('Determine principal component time series %d',k))
        for i=0:n-m
            %modify for nan
            prod=x(i+1:i+m).*ro(:,k);
            igood=find(~isnan(prod));
            ngood=length(igood);
            %must have at least m*f good points
            if ngood>=m*f 
                xk(k,i+1)=sum(prod(igood))*m/ngood;
                xk_(k,i+1)=sum(prod(igood)); %rs, no scaling here
            end
        end
    end
end

%--------------------------------------------------------------------------

function xf=rc(xk,ro)
%function xf=rc(xk,ro)
%
% purpose: make normalized time series (reconstructed components)
%          with vectorized commands for more speed
%
% inputs: xk=principal components 
%         ro=eigenvectors
%
% outputs: xf=reconstructed components (normalized)
%
% by David Schoellhamer 9/24/93
%
%
m=max(size(ro));
n=max(size(xk))+m-1;
nk=min(size(ro));
xf=zeros(n,nk);
%
% use vectorized Vautard and Ghil eqn B.5
%
%for ik=1:nk
  for i=1:n
    if i==1
      xf(i,:)=(rot90(xk(:,1:i)).*ro(1:i,:))/i;
    elseif i<m
      xf(i,:)=sum(rot90(xk(:,1:i)).*ro(1:i,:))/i;
    elseif i==n
      xf(i,:)=(rot90(xk(:,i-m+1:n-m+1)).*ro(i-n+m:m,:))/(n-i+1);
    else if i>n-m+1
      xf(i,:)=sum(rot90(xk(:,i-m+1:n-m+1)).*ro(i-n+m:m,:))/(n-i+1);
    else
      i1=i-m+1;
      xf(i,:)=sum(rot90(xk(:,i1:i)).*ro(:,:))/m;
    end
%    if rem(i,1000)==0 %rs commented out
%      disp(sprintf('%d of %d data points reconstructed',i,n))
%    end
  end
%end
end
xf=xf';
end

%--------------------------------------------------------------------------

function x_in_cont=getdata()

%read a trend and jumps removed gappy hector file from Hadi and return with NaNs in gaps
%0=ew,1=NS,2=ud, ordered but not uniformly spaced, column vector
filename_root='~/CentralUS_Ubuntu_toSSA/Stns_Dir/';
sitename='1LSU';
compname='E';
compnum='0';
filename='1LSU_E/1LSU_0_SSA.dat'; %(all real) input data has gaps
full_filename=strcat(filename_root,sitename,"_",compname,"/",sitename,"_",compnum,"_SSA.dat");

input_data_discont=load(full_filename);

x_in_discont=input_data_discont(:,3); %discontinuous in time data vector, cols: fracyr, mjd, position (1 comp/file)
x_in_discont_indx=input_data_discont(:,1)-input_data_discont(1,1)+1; %this has indices of days with data
N=input_data_discont(end,1)-input_data_discont(1,1)+1; %num days from first to last, >= to number days with data
x_in_cont=nan(N,1); %to make a copy of input vector that is continuous, fill with NaNs, will remain in holes
x_in_cont(x_in_discont_indx)=x_in_discont; %move in data, discont_in_indx has index
end

