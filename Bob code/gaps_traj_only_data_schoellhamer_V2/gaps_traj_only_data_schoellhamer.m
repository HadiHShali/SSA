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
x_mean=mean(x_in);
x_std=std(x_in);
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
m=fix(n_data/2); %can only have as many RC as [fix(number data/2)]. m≤half number data (and continuous length≥ndata)
%this is (generally) max embedding dimension (max number RCs and shortest
%sections of time series)
%m=1278; %set embedding dimension to 3.5 years (rounded to integer), gives number of columns in trajectory matrix and number RCs we will obtain
l=n_cont_time-m+1; %number rows in trajectory matrix (for m=n/2, with n even, is one taller than wide, with n odd is 2 taller than wide), as m is smaller, l is longer

colind=1:m; %colmn index vector
rowind=0:l-1; %row index vector
trajmatind=colind+rowind'; %vectorized build of trajectory matris
trajmat=x(trajmatind);

OutStrucNaNs=SSABasicCleanUCLA(trajmat); %locally defined below - fixed for missing data

%how to pick this automatically
% how to pick this automatically????
maxeig=14;
%maxeig=18;
%maxeig=29;
maxeig=50;
%maxeig=100;
%maxeig=m/2;
%maxeig=2;
go=1;

while go %0 is false, anything else true

    figname=strcat("eigenvalues, max eig=",num2str(maxeig));
    figno=newfig(figno,figname);
    semilogy(OutStrucNaNs.LAMBDA,'b+-')
    hold on
    semilogy(OutStrucNaNs.LAMBDA(1:maxeig),'ro-')
    closefig(figname)
    
    RC=OutStrucNaNs.RCunscaled;
    %RC=OutStrucNaNs.RCscaled;
    RCNaNs=InsertNaNs(RC,NaNindx);
    reconsNaNs=sum(RCNaNs,2);
    diff_ck(x,reconsNaNs,"x-reconsNaNs")
    %x_reconsNaNs=reconsNaNs*x_std_omitnan+x_mean_omitnan;
    %diff_ck(x_in,x_reconsNaNs,"x_in-x_reconsNaNs")
    
    figname=strcat("recon gappy: Ntime=", num2str(n_cont_time)," Ndata=",num2str(n_data),...
        ", M=", num2str(m), ", L=", num2str(l), ", max eig=", num2str(maxeig));
    figno=newfig(figno,figname);
    plot(t,x,'b')
    hold on
    plot(t,reconsNaNs,'r--')
    rc=sum(RCNaNs(:,[1:maxeig]),2);
    rc(NaNindx)=NaN;
    plot(t,rc,'k')
    xlabel('time (days)')
    ylabel('position')
    legend({'in' 'full RC' 'gappy sig model'},'Location','southeast')
    %closefig(figname)
    
    figname=strcat("recon gappy noise TS, max eig=", num2str(maxeig));
    figno=newfig(figno,figname);
    ruidoest=sum(RCNaNs(:,[maxeig+1:end]),2);
    plot(t,ruidoest,'r')
    closefig(figname)
    
    %go=input(strcat("current maxeig=", num2str(maxeig)," enter new maxeig, <CR> for quit "));
    go=0;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function OutStruc = SSABasicCleanUCLA(trajmat)
    disp('entering local SSABasicCleanUCLA')
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
 
    Csave=C;

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

        %get number "good" terms in trajmat*RHO, all terms in RHO are good, so only need to know number bad terms in trajmat
        numterms=sum(trajmat_ele_TF,2); %number good terms in each row, max for good row is number columns
        scl=M./numterms; %this must be >= to 1;
        minscl=min(scl(:)); %probably don't need to do this check
        if minscl<1
            disp('problem with scale in making PCs')
            return
        end        
        PCunscaled=trajmat_ele_good*RHO;
        PCscaled=(scl.*trajmat_ele_good)*RHO; %Nx1 vector times NxM matrix with singleton expansion
        disp('scaled and unscaled, has NaNs') %can't figure out this part, unscaled seems to work
        diff_ck(PCunscaled,PCscaled,'PCunscaled-PCscaled, should be different')
    else
        PC = trajmat_ele_good*RHO; %make principal components, will be size of trajectory matrix as RHO is square
        disp('UNSCALED, no NaNs');
    end


    % Calculate reconstructed components RC
    % In order to determine the reconstructed components RC, we have to invert
    % the projecting PC = Y*RHO; i.e. RC = Y*RHO*RHO'=PC*RHO'. Averaging along
    % anti-diagonals gives the RCs for the original input X.

    tstscaled=sum(isnan(PCscaled(:)));
    tstunscaled=sum(isnan(PCunscaled(:))); %don't have to test PC
    if tstscaled~=0 | tstunscaled~=0
        disp(strcat("missing pc terms - have to fix, PCunscaled ",num2str(tstunscaled),...
            " PCscaled ",num2str(tstscaled) ))
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

    RCscaled=zeros(N,M);
    RCunscaled=RCscaled;
    tic
    for m=1:M
        bufscaled=PCscaled(:,m)*RHO(:,m)'; % invert projection
        bufscaled=flipud(bufscaled);
        bufunscaled=PCunscaled(:,m)*RHO(:,m)'; % invert projection
        bufunscaled=flipud(bufunscaled);

        for n=1:N % anti-diagonal averaging
                %here the rescaling gets done by having NaNs in the
                %diagonal sum - that provides a NaN answer for the (n,m) element
                %of the RC
                RCunscaled(n,m)=mean( diag(bufunscaled,-(N-M+1)+n)); %this takes forever but I've not been able to beat it (yet)
                RCscaled(n,m)=mean( diag(bufscaled,-(N-M+1)+n)); %this takes forever but I've not been able to beat it (yet)
        end
    end
    toc

    diff_ck(RCunscaled,RCscaled,'RCunscaled-RCscaled')

    OutStruc.RCscaled=RCscaled;
    OutStruc.RCunscaled=RCunscaled;
    OutStruc.LAMBDA=LAMBDA;
    OutStruc.RHO=RHO;
    OutStruc.num_nans=num_nans;
    disp('leaving local SSABasicCleanUCLA')

end

