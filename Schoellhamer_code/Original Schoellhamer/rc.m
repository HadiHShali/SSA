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
    if rem(i,1000)==0
      disp(sprintf('%d of %d data points reconstructed',i,n))
    end
  end
%end
end
xf=xf';