%% M-SSA varimax tutorial
% This Matlab(R) tutorial demonstrates the application of the varimax algorithm
% to the eigenvectors of a multichannel singular spectrum analysis (M-SSA).
%
%   Requires varimax.m
%
% For more mathematical details see
%
% # Groth, A. & Ghil, M. (2011), 'Multivariate singular spectrum analysis
% and the road to phase synchronization', Physical Review E 84, 036206.
%
% For more tutorials see <https://dept.atmos.ucla.edu/tcd/matlab-tutorials>
%
% _Copyright (c) 2010-2016, Andreas Groth, University of California, Los
% Angeles. All rights reserved. Redistribution and use in source and binary
% forms, with or without modification, are permitted provided that the
% following conditions are met: (1) Redistributions of source code must
% retain the above copyright notice, this list of conditions and the
% following disclaimer. (2) Redistributions in binary form must reproduce
% the above copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials provided with the
% distribution._

%% Generate time series
% Generate |D=6| time series of length |N=300|, each of them a linear combination of
% four sinusoids with different frequency |f0| and variance |f0Var|.

N=300;
t=(1:N)';

% fundamental frequencies
f0 = 1./[6 8 10 13];

% variance
f0Var = [
         0.2  0.0 0.3  0.4 ;...         
         0.4  0.3 0.4  0.0 ;...
         0.2  0.3 0.0  0.4 ;...
         0.0  0.4 0.4  0.2 ;...
         0.2  0.4 0.0  0.4 ;...
         0.3  0.0 0.5  0.4];
D=size(f0Var,1);

% combination of sinusoids
x=zeros(N,D);
for d=1:D
  for pos=1:length(f0)
    x(:,d)=x(:,d) + sqrt(f0Var(d,pos)) * sin( 2*pi*f0(pos)*t+rand(1)*2*pi );
  end
end

%% M-SSA analysis
% Perform M-SSA with the Broomhead-King approach of a covariance-matrix
% estimation |C|. The M-SSA window length is |M=40|. 

M=40;

% time-delayed embedding
xtde=zeros(N,N,D);
for d=1:D
  xtde(:,:,d)=hankel(x(:,d));
end
xtde=xtde(1:N-M+1,1:M,:);
xtde=reshape(xtde,N-M+1,D*M,1);

% M-SSA analysis
C=xtde'*xtde/(N-M+1); % Broomhead and King (1986)
[EV,EW]=eig(C);
EW=diag(EW);
[EW ind]=sort(EW,'descend');
EV=EV(:,ind);

%%
% Plot the leading eight eigenvectors arranged as columns in |EV|. The
% eigenvectors have length |D*M|.
figure(1)
for pos=1:2*length(f0)
  subplot(2*length(f0),1,pos)
  plot(reshape(1:D*M,M,D),reshape(EV(:,pos),M,D))
  xlim([0 D*M]);
  ylim([-1 1]*0.3);
  set(gca,'XTick',0:M:D*M,'XGrid','on');
  ylabel(sprintf('{\\bf e}_%d',pos),'Rotation',0);
end
xlabel('Embedding dimension');
subplot(2*length(f0),1,1);
title('Unrotated eigenvectors');

%% Varimax rotation
% Next, we apply varimax rotation to the first |S=20| eigenvectors. The
% eigenvectors |EV| are scaled by their singular value |sqrt(EW)| and then
% reshaped. With this reshaping, the function |varimax.m| can distinguish PCA
% eigenvectors from M-SSA eigenvectors. 

S=20; % number of rotated eigenvectors

EVscaled=EV(:,1:S)*diag(sqrt(EW(1:S))); % scaling
EVreshape=reshape(EVscaled,M,D,S);      % reshaping is necessary for varimax.m

[dummy T]=varimax(EVreshape);           % varimax rotation

EV(:,1:S)=EV(:,1:S)*T;                  % update eigenvectors
EW(1:S)=diag(T'*diag(EW(1:S))*T);       % update eigenvalues

% sort again eigenelements
[EW ind]=sort(EW,'descend');
EV=EV(:,ind);

%%
% Plot the leading eight eigenvectors in |EV| after varimax rotation.
figure(2)
for pos=1:2*length(f0)
  subplot(2*length(f0),1,pos)
  plot(reshape(1:D*M,M,D),reshape(EV(:,pos),M,D))
  xlim([0 D*M]);
  ylim([-1 1]*0.3);
  set(gca,'XTick',0:M:D*M,'XGrid','on');
  ylabel(sprintf('{\\bf e}_%d',pos),'Rotation',0);
end
xlabel('Embedding dimension');
subplot(2*length(f0),1,1);
title('Rotated eigenvectors');