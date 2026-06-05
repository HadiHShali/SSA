function OutStruc = SSABasicCleanUCLAtoeplitz(varargin)
%function OutStruc = SSABasicCleanUCLAtoeplitz(trajmat,Ctoep)
    disp(strcat("entering local SSABasicCleanUCLAtoeplitz, nargin ",num2str(nargin)))
    trajmat=cell2mat(varargin(1));
    if nargin==2
        toep=cell2mat(varargin(2));
        disp("use toeplitz for covariance")
    else
        disp("use traj'*traj for covaraince")
    end
    %disp(InStruc.RCs)
    %from%https://www.mathworks.com/matlabcentral/fileexchange/58967-singular-spectrum-analysis-beginners-guide
    
    %input structure InStruc.FIELD fields - UON are required
    %       trajmat - the trajectory matrix, is embedding dimension wide
    %       and L tall where the total number of points in the input vector
    %       is N and N=L+M-1. It has been pre-prepared witn NaN elements
    %       for missing data
    %       that will be replaced with 0 (actual data will never be equal to floating point
    %       zero)
    %       num_nans is the number of NaNs in the inut vector

    [L,M]=size(trajmat); %M is embedding dimension (shold be based on number of data, i.e. excluding NaNs),
    %L is length of the segments and these segments include NaNs in the
    %count
    
    N=L+M-1; %N is total number points in data vector (spans time with daily sampling from 1st to last day) including NaNs,
    %M is embedding dimension (data considered "gap" free?) - some confusion
    %about embedding dimension - is it on full time series or only days
    %with data???? I'm using full time series. If using something other than half of length of time
    %series less number NaNs - what does it represent, actual time? what is requirement wrt number of
    %days with data and total number of days (is max data/2 or total
    %days/2?).

    originalvec=[trajmat(1:L,1); trajmat(end,2:end)'];
    num_nans=sum(isnan(originalvec));
    num_data=N-num_nans;

    if num_nans~=0 %have to compute trajmat'*trajmat using 0s for missing data (NaNs) and adjust for number of data (non NaNs) in dot products
        disp('hay NaNs');
        hay_nans=1;
        trajmat_ele_good=trajmat; %put 0's into elements with NaN (this is a copy to work on, don't clobber original)
        trajmat_ele_good(isnan(trajmat_ele_good))=0; %all NaNs replaced with zeros
        trajmat_ele_TF=ones(size(trajmat)); %traj mat element good, 1=true (has data), 0=false (O for Nan)
        trajmat_ele_TF(isnan(trajmat))=0;
        pct_good=sum(trajmat_ele_TF(:))/length(trajmat(:)) %size(trajmat(:)) returns a 1x2 vector of length(trajmat(:)) by 1. could use size(trajmat(:),1),
        %from command line size(trajmat(:)) returns single value, but [a b]=size(trajmat(:)) returns a=number points and b=1

        [indxelegood,logical]=find(~isnan(trajmat)); %linear index for good data
        [indxeleNaN,logical]=find(isnan(trajmat)); %linear index for missing data, these 2 sets of indices should not overlap and their union should cover everything

        %it is hard to get 0 elements in the following two products for
        %resonable time series.
        %remembering that the elements of matrix multiplication are the dot
        %products of the rows and columns as vectors the only way to get 0 is
        %for all elements to be zero or the vectors being perpendicular. 
        %NaNadj gives the number of non-zero products in the
        %dot product so zero elements in NaNadj signify that the vectors are either both 0
        %or always opposite (when an element is 1 in the row, the corresponding one in the column is 0).
        %Since the first and last element of the time series cannot be NaN=0,
        %there will always be at least 2 elements of NaNajd that are not zero.
        %see sets of pathological series at the beginning - probably will
        %not arise, but check!

        %Note - to pull elements out of matrix using vectors to specify indices
        %have to convert to linear indices. A(sub2ind(size(A),r,c)), where r has
        %the 1st index and c has the 2nd index.

        %For C one can get zero if all the products are zero or the vectors are perpendicular.
        %so some 0 elements in C can be "valid" (but they are kept as the
        %NaNadj will not catch them so the M/n_used will be correct.

        C=trajmat_ele_good'*trajmat_ele_good; %this creates the covariance matrix from the trajectory matrix, missing points have value 0 and don't contribute, have to "normalize" by number good points in each dot product
        %in general - elements of C are not zero (need whole column/row=0 or perpendicular vectors (have to try with monster gaps)
        NaNadj=trajmat_ele_TF'*trajmat_ele_TF; %trajmat_ele_TF has true (1) for data, and false (0) for missiong data, this matrix product gives number good points in each dot product - for normalization
        %(need whole column/row=0 or perpendicular vectors (have to try with monster gaps)
        C=C./NaNadj; %Normalize by number terms in dot product. If any NaNadj are zero the result will be Inf, check - if this happens what do we do??? wait till it happens and figure it out (this should only happen if any data point is exactly zero).

        problem=sum(isinf(C(:))); %check for problem, advise, and stop
        if problem>0
            disp(strcat(num2str(problem)," bad elements in C (infinity, divide by 0 correlation)"))
            return
        end
    else
        hay_nans=0; %if no nans do normally
        C=trajmat'*trajmat;
    end
 
    if nargin==2
        C=toep;
        disp('use toeplitz')
    end
    Csave=C;

%here we have the lagged autocorrelation matrix in C (equivalent to eq 4 in
%Schoellmamer)


    % Calculate eigenvalues LAMBDA and eigenvectors RHO
    % In order to determine the eigenvalues and eigenvectors of C, we use the
    % function EIG. This function returns two matrices, the matrix RHO with
    % eigenvectors arranged in columns, and the matrix LAMBDA with eigenvalues
    % along the diagonal.

    [RHO,LAMBDAM] = eig(C);
    RHOexist=sum(isnan(RHO));
    if RHOexist>0
        disp('no eigenvalue/eigenvector solution')
        return
    end
    LAMBDA = diag(LAMBDAM);               % extract the diagonal elements
    LAMBDA(LAMBDA<0)=1e-16; %get rid of illegal negative eigenvalues due to numerical instabilities (very rare)
    [LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
    RHO = RHO(:,ind);                    % and eigenvectors (using indices from sort)

    %the eigenvectors do not contain missing data, so we only have to
    %adjust for what they are multiplied by

    % Calculate principal components PC
    % The principal components are given as the scalar product between Y, the
    % time-delayed embedding of X, and the eigenvectors RHO taking into 
    % account the missing data NaNs by replacing them with zeros and then
    % adjusting the ave of the row column dot product by the number of good terms.
   

    if hay_nans %have to modify PC = trajmat*RHO for missing data [NaNs (=0)]
        disp('hay NaNs');
    end

    PC = trajmat_ele_good*RHO; %make principal components, will be size of trajectory matrix as RHO is square

    % Calculate reconstructed components RC
    % In order to determine the reconstructed components RC, we have to invert
    % the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
    % anti-diagonals gives the RCs for the original input X.

    tst=sum(isnan(PC(:)));
    if tst~=0
        disp(strcat("missing pc terms - have to fix, PC ",num2str(tst)))
        return
    %if any PC component (element in PC, not complete PC vector?) is missing, then the RC value will be missing
    %this is confusing - I think it means one element of the RC being
    %calculated is missing. PC==0 only when there is a set of 0s (were NaNs in time series) in the
    %trajectory matrix of length>M, then the whole ROW of the PC matrix is 0.
    %I think we change the 0's to NaNs so in the next matrix multiply those dot products will be NaN
    %then do the diagonal averaging on buf, omitting those NaNs.
    %this keeps the RC the correct size and weighted correctly and things still line up at the end
    %worry about it when it happens
    end
    disp('no NaNs in PCs'); %in either case (gaps or no gaps) there should not be NaNs in the PCs

    tic
    for m=1:M

        buf=PC(:,m)*RHO(:,m)'; % invert projection
        buf=flipud(buf);

        for n=1:N % anti-diagonal averaging
                %here the rescaling gets done by having NaNs in the
                %diagonal sum - that provides a NaN answer for the (n,m) element
                %of the RC
                RC(n,m)=mean( diag(buf,-(N-M+1)+n)); %this takes forever but I've not been able to beat it (yet)
        end
    end
    toc

    OutStruc.RC=RC;
    OutStruc.LAMBDA=LAMBDA;
    OutStruc.RHO=RHO;
    OutStruc.num_nans=num_nans;
    disp('leaving local SSABasicCleanUCLA')

end

