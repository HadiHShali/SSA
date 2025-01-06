function [xmean,xstd,xk,ro,lam]=vssanan(ts,m,f)
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
disp('Compute autocorrelation function')
for j1=1:m
%  summ=0;
  j=j1-1;
%  for i=1:n-j
%    summ=summ+x(i)*x(i+j);
%  end
%  c(j1)=summ/(n-j-1);
%
% modify for nan
  prod=x(1:n-j).*x(j+1:n);
  igood=find(~isnan(prod));
  c(j1)=sum(prod(igood))/(length(igood)-1);
end
%
% form the covariance matrix a, divide by m (eqn b.3)
%
disp('Form covariance matrix a')
%a=diag(c(1)*ones(m,1));
%for j=2:m
%  a=a+diag(c(j)*ones(m+1-j,1),j-1)+diag(c(j)*ones(m+1-j,1),1-j);
%end
%a=a/m;
a=toeplitz(c)/m;
%
% determine eigenvalues and eigenvectors of a (solve b.3)
%
disp('Determine eigenvalues and eigenvectors of a');
[z,eval]=eig(a);
lam=diag(eval);
%
% sort by descending eigenvalues
%
disp('Sort by descending eigenvalues');
[lam,ilam]=sort(lam);
lam=flipud(lam);
ilam=flipud(ilam);
summ=0;
kmax=0;
semilogy(lam(1:min([m 30])),'o')
nkmax=0;
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
  disp(sprintf('Eigenvalue of pc %d is %g, sum of variance is %g',...
               kmax,lam(kmax),summ))
end
%
% determine principal component time series (eqn b.6)
%
xk=nan*ones(kmax,n-m+1);

for k=1:kmax
  disp(sprintf('Determine principal component time series %d',k))
  for i=0:n-m
%   modify for nan
    prod=x(i+1:i+m).*ro(:,k);
    igood=find(~isnan(prod));
    ngood=length(igood);
%   must have at least m*f good points
    if ngood>=m*f 
      xk(k,i+1)=sum(prod(igood))*m/ngood;
    end
  end
end
