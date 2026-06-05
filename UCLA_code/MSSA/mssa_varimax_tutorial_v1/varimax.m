function [B T crit] = varimax(A,reltol,maxit,normalize,G)
%VARIMAX Varimax rotation of eigenvectors
%   B=VARIMAX(V) returns the result of an orthogonal varimax rotation of
%   the column vectors in A.
%
%   [B,T,CRIT]=VARIMAX(...) also returns the rotation matrix B=V*T and in
%   CRIT the variance at each iteration step.
%
%   Additional options:
%
%   [..]=VARIMAX(V,RELTOL,MAXIT) sets the break condition for CRIT<RELTOL,
%   and the maximal number MAXIT of iterations. (default values
%   RELTOL=sqrt(eps(class(V))) , MAXIT=1000)
%
%   [..]=VARIMAX(V,RELTOL,MAXIT,NORM) switches normalization of rows in V
%   on/off; i.e. NORM=true/false (default NORM=true).
%
%   [..]=VARIMAX(V,RELTOL,MAXIT,NORM,G) additionally restricts the
%   calulation of the varimax criterion between groups, defined in G. The
%   numeric vector G of length D, defines groups in V between which the
%   varimax criterion is maximized. Groups with negative values in G are
%   excluded from the criterion. Leave empty for no grouping (default
%   value=[]).
%
%   Especially for the application on M-SSA eigenvectors, V can have a third
%   dimension. In this case, each space-time eigenvector is a matrix along
%   the first two dimensions (time along first dimension and space along
%   second dimension), and the eigenvectors are arranged along the third
%   dimension. Then, the algorithm ensures that the variance is not
%   maximized within time but only within space.
%
%   Note that reshaping is a special case of grouping with G. The
%   grouping requires no reshaping of the M-SSA eigenvectors, with
%     G = kron((1:D),ones(1,M))
%
%   Example for M-SSA eigenvectors:
%     % Let's assume that we obtain the eigenvectors V from M-SSA such as
%     D=10; M=5;
%     X=randn(100,D);
%     [RC,A,V]=mssa(X,M);
%     % Then, we have to reshape V before passing it to VARIMAX
%     V2=reshape(V(:,1:S),M,D,S);
%     B=VARIMAX(V2);
%     % Note that we rotate only the first S eigenvectors in this example
%
% See also MSSA, ROTATEFACTORS.

% Andreas Groth, ENS, Paris, France and UCLA, Los Angeles.
% 07/08/2010 First implemenation
% 03/23/2016 Added grouping option
% 03/31/2016 Added normalization option

% Parts are taken from Matlab function rotatefactors.m

if nargin < 2
  reltol=[];
end
if nargin < 3
  maxit=[];
end
if nargin < 4
  normalize=true;
end
if nargin < 5
  G =[];
end

B=A;

if length(size(B))>2 % ------ M-SSA --------
  [M,D,S]=size(B);
  % normalize along parts of eigenvectors
  if normalize
    h=sum(sum(B.^2,1),3);
    if ~isempty(G) % replace h by sum of h in each group of G
      Gs=bsxfun(@eq,G(:),G(:)');
      h=h*Gs;
    end
    B=bsxfun(@times,B,1./sqrt(h));
  end
  B=reshape(B,D*M,S,1);
  D=D*M;
else       % -------- Classical PCA --------
  [D,S] = size(B);
  M=1;
  % normalize along each row
  if normalize
    h=sum(B.^2,2);
    if ~isempty(G) % replace h by sum of h in each group of G
      Gs=bsxfun(@eq,G(:),G(:)');
      h=Gs*h;
    end
    B=bsxfun(@times,B,1./sqrt(h));
  end
end;

if isempty(reltol)
  reltol = sqrt(eps(class(B)));
end
if isempty(maxit)
  maxit = 1000;
end
if ~isempty(G)
  Gc=unique(G);
  Gc(Gc<0)=[];
  G=bsxfun(@eq,G(:),Gc(:)');
end

T=eye(S);

if nargout>2
  crit=nan(maxit,1);
end;

% Use a sequence of bivariate rotations
for iter = 1:maxit
  maxTheta = 0;
  for i = 1:(S-1)
    for j = (i+1):S
      u = B(:,i).^2 - B(:,j).^2;
      v = 2*B(:,i).*B(:,j);
      % Consider for each eigenvector-part a single value = sum
      if M>1
        u= sum(reshape(u,M,[]))';
        v= sum(reshape(v,M,[]))';
      end;
      % only enhance differences between groups
      if ~isempty(G)
        u=G'*u;
        v=G'*v;
      end
      
      usum = sum(u,1);
      vsum = sum(v,1);
      numer = 2*u'*v - 2*usum*vsum/D;
      denom = u'*u - v'*v - (usum^2 - vsum^2)/D;
      theta = atan2(numer, denom) / 4;
      maxTheta = max(maxTheta, abs(theta));
      Tij = [cos(theta) -sin(theta); sin(theta) cos(theta)];
      B(:,[i,j]) = B(:,[i,j]) * Tij;
      T(:,[i,j]) = T(:,[i,j]) * Tij;
    end
  end
  if nargout>2
    if M>1
      dummy=squeeze(sum(reshape(B,M,[],S).^2));
    else
      dummy=B.^2;
    end;
    crit(iter)=sum(var(dummy,1));
  end;
  if (maxTheta < reltol)
    break;
  end
end
if nargout<3 && iter==maxit
  warning('Convergence limit not reached!');
end

% unnormalize
if normalize
  if M>1
    B=reshape(B,M,D/M,S);
    B=B.*repmat(h,[M 1 S]);
    B=reshape(B,D,S,1);
  else
    B=B.*repmat(h,1,S);
  end;
end