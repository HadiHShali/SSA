%from https://dept.atmos.ucla.edu/tcd/ssa-tutorial-matlab
%--------------------------------------------------------------------------
%...
%--------------------------------------------------------------------------
clear
close all

%% Set general Parameters
M = 30;    % window length = embedding dimension
N = 200;   % length of generated time series
T = 22;    % period length of sine function
stdnoise = 1; % noise-to-signal ratio


%% Create time series X
%First of all, we generate a time series, a sine function of length N with observational white noise

t = (1:N)';
Xsig = sin(2*pi*t/T);
noise = stdnoise*randn(size(Xsig));
X = Xsig + noise;
X = X - mean(X);            % remove mean value
X = X/std(X,1);             % normalize to standard deviation 1

%X=1:10;
%NN=length(X);
NN=N;
X=X(1:NN);
%MM=fix(NN/2);
MM=M;
figure;
set(gcf,'name','Time series X');
%clf;
plot(t(1:NN),X,'b-');
hold on
plot(t(1:NN),Xsig,'r-');
legend('Signal+noise','Original Signal')
title('Time sereis')


%% Calculate covariance matrix C (Toeplitz approach)
%Next, we calculate the covariance matrix. There are several numerical approaches to estimate C. Here, we calculate the covariance
% function with CORR and build C with the function TOEPLITZ.

covX = xcorr(X,M-1,'unbiased'); %none for different lengths, Rxy(m) biased 1/N, Rxy(m) unbiased 1/(N-|m|), normalized - max 1)
Ctoep=toeplitz(covX(M:end));

%--------------------------------------------------------------------------

%test vectorized (above) against "real" example of format/matfor - and using caps
%only this does not interfere with sqrt -1!! all papers have this abomination
ckpaper=1;
ckpaper=0;
if ckpaper
    for J=0:MM-1 %do 
    CS(J+1)=0;
    for I=1:NN-J
        CS(J+1)=CS(J+1)+X(I)*X(I+J);
    end
    CS(J+1)=CS(J+1)/(NN-J);
end
CS=CS';
max(abs(CS-Ctoep(:,1))) %only need to check against 1st column of toeplitz matrix
%CTOEP=toeplitz(CS); %make toeplitz matrix (matlab vectorized), could do
%with complicated format/matfor loops
end
%--------------------------------------------------------------------------

figure;
set(gcf,'name','Covariance matrix');
clf;
imagesc(Ctoep); %by default imagesc inverts the y axis, positive down. This makes
%matrices, and things like landsat images, start at upper left and go across and down
%x(1 1) x(1 2)
%x(2 1) x(2 2)
%as in Math notation for matrices, not as x axis positive right and y axis
%positive "up" page for plotting
axis square
%set(gca,'clim',[-1 1]);
colorbar
title('Covariance by Toeplitz approach')

%% Calculate covariance matrix (now by trajectory approach)
%An alternative approach is to determine C directly from the scalar product of Y, the time-delayed embedding of X. Although this estimation of C does not give a Toeplitz structure, with the eigenvectors not being symmetric or antisymmetric,
% it ensures a positive semi-definite covariance matrix.
N = length(X);
Y=zeros(N-M+1,M);
for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end

%--------------------------------------------------------------------------
%above can also be vectorized
%--------------------------------------------------------------------------

Cemb=Y'*Y / (N-M+1);
figure;
set(gcf,'name','Covariance matrix');
clf;
imagesc(Cemb);
axis square
set(gca,'clim',[-1 1]);
colorbar
title('Covariance by trajectory approach')

%% Choose covariance estimation
%Choose between Toeplitz approach (cf. Vautard & Ghil) and trajectory approach (cf. Broomhead & King).
% C=Ctoep;
C=Cemb;

%% Calculate eigenvalues (LAMBDAs) and eigenvectors (RHOs)
%In order to determine the eigenvalues and eigenvectors of C, we use the function EIG. This function returns two matrices,
% the matrix RHO with eigenvectors arranged in columns, and the matrix LAMBDA with eigenvalues along the diagonal.

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors
figure;
set(gcf,'name','Eigenvectors RHO and Eigenvalues LAMBDA')
clf;
subplot(3,1,1);
plot(LAMBDA,'o-');
title('Eigenvalues')
subplot(3,1,2);
plot(RHO(:,1:2), '-');
legend('1', '2');
title('Eigenvectors')
subplot(3,1,3);
plot(RHO(:,3:4), '-');
legend('3', '4');

%% Calculate principal components PC
%The principal components are given as the scalar product between Y, the time-delayed embedding of X, and the eigenvectors RHO.

PC = Y*RHO;
figure;
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t(1:N-M+1),PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end
sgtitle('Principal Components (PCs)');

%% Calculate reconstructed components RC
%In order to determine the reconstructed components RC, we have to invert the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'.
% Averaging along anti-diagonals gives the RCs for the original input X.
RC=zeros(N,M);
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end
figure;
set(gcf,'name','Reconstructed components (RCs)')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
end
sgtitle('Reconstructed Components (RCs)')

%% Compare reconstruction and original time series
%Note that the original time series X can be completely reconstructed by the sum of all reconstructed components RC (upper panel).
% The sine function can be reconstructed with the first pair of RCs (lower panel).

figure;
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(2,1,1)
plot(t,X,'b',t,sum(RC(:,:),2),'r-.');
legend('Original','Complete reconstruction');
title('Compare reconstruction and original time series')
subplot(2,1,2)
plot(t,X,'b','LineWidth',2);
plot(t,X,'b-',t,sum(RC(:,1:2),2),'r-');
legend('Original','Reconstruction with RCs 1-2');
