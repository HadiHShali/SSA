# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:21:31 2025

@author: GeodesyLab
"""
# %matplotlib qt

# %% Clear workspace and setup
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz  # For Toeplitz matrix generation
import os
import pandas as pd

# SSA function
# % SSA Base Function
def SSABasicCleanUCLAtoeplitz(trajmat, toep=None, maxieg=None):
    """
    Perform Singular Spectrum Analysis (SSA) to clean the trajectory matrix with or without Toeplitz covariance.

    Parameters:
    - trajmat: Trajectory matrix (embedding dimension wide, segments tall).
    - toep: Optional Toeplitz matrix for covariance, if not provided, will use trajmat' * trajmat.

    Returns:
    - OutStruc: Dictionary with reconstructed components (RC), eigenvalues (LAMBDA), eigenvectors (RHO), and number of NaNs.
    """
    print(f"Entering local SSABasicCleanUCLAtoeplitz, number of inputs: {len(locals())}")

    # Check the input dimensions and prepare the trajectory matrix
    L, M = trajmat.shape  # L is the number of segments, M is the embedding dimension

    # Handle missing data (NaNs) in the trajectory matrix
    originalvec = np.concatenate((trajmat[:, 0], trajmat[-1, 1:]))  # Concatenate first column and last row (excluding first element)
    num_nans = np.sum(np.isnan(originalvec))  # Count the number of NaNs in the original vector
    N = L + M - 1  # Total number of data points (L + M - 1)
    num_data = N - num_nans  # Number of data points excluding NaNs
    
    if num_nans != 0:
        print("NaNs detected in the input matrix.")
        trajmat_ele_good = np.copy(trajmat)
        trajmat_ele_good[np.isnan(trajmat_ele_good)] = 0  # Replace NaNs with 0s
        trajmat_ele_TF = np.ones_like(trajmat)  # Create a true/false mask for data points
        trajmat_ele_TF[np.isnan(trajmat)] = 0
        
        # Percentage of good data
        pct_good = np.sum(trajmat_ele_TF) / trajmat.size
        print(f"Percentage of good data: {pct_good:.4f}")
        
        # Indices for non-NaN and NaN elements
        indxelegood = np.where(~np.isnan(trajmat))  # Indices for non-NaN data
        indxeleNaN = np.where(np.isnan(trajmat))  # Indices for NaN data

        # Covariance matrix using good data (ignoring NaNs)
        #%this creates the covariance matrix from the trajectory matrix, missing points have value 0 and don't contribute, have to "normalize" by number good points in each dot product
        #in general - elements of C are not zero (need whole column/row=0 or perpendicular vectors (have to try with monster gaps)
        C = (trajmat_ele_good.T) @ trajmat_ele_good
        
        #trajmat_ele_TF has true (1) for data, and false (0) for missiong data, this matrix product gives number good points in each dot product - for normalization
        NaNadj = np.dot(trajmat_ele_TF.T, trajmat_ele_TF)  # gives the number of non-zero products in the dot product. 

        C = C/NaNadj  # Normalize by number of good terms

        # Check for problems in the covariance matrix (inf values)
        if np.any(np.isinf(C)):
            print("Bad elements in covariance matrix (infinity or divide by zero)")
            return

    else:
        print("No NaNs detected, using normal covariance")
        C = np.dot(trajmat.T, trajmat)  # Standard covariance calculation if no NaNs
    
    # If Toeplitz matrix is provided, use it instead of calculated covariance
    if toep is not None:
        C = toep
        print("Using Toeplitz matrix for covariance")
    
    # Calculate eigenvalues (LAMBDA) and eigenvectors (RHO)
    try:
        LAMBDAM, RHO = np.linalg.eig(C)  # Eigen decomposition
    except np.linalg.LinAlgError:
        print("Eigenvalue/eigenvector calculation failed.")
        return
    
    # Handle negative eigenvalues due to numerical issues
    LAMBDA = LAMBDAM
    #LAMBDA = np.copy(LAMBDA)
    LAMBDA[LAMBDA < 0] = 1e-16  # Set negative eigenvalues to a small positive number
    
    sorted_indices = np.argsort(LAMBDA)[::-1]  # Get indices for sorting in descending order

    [LAMBDA, ind] =[ np.sort(LAMBDA)[::-1], np.argsort(LAMBDA)[::-1]]  # Sort eigenvalues in descending order
    RHO = RHO[:, ind]  # Reorder eigenvectors according to sorted eigenvalues

    # Calculate principal components (PC)
    #    % The principal components are given as the scalar product between Y, the
    # time-delayed embedding of X, and the eigenvectors RHO taking into 
    # account the missing data NaNs by replacing them with zeros and then
    # adjusting the ave of the row column dot product by the number of good terms.
    # if num_nans > 0:
    #     print("NaNs detected in principal components, need to handle them.")
    

    if maxieg is not None:
        maxieg = min(maxieg, len(LAMBDA))
        LAMBDA_maxieg = LAMBDA[:maxieg]
        RHO_maxieg = RHO[:, :maxieg]
    
        PC = np.dot(trajmat_ele_good, RHO_maxieg)
        if np.any(np.isnan(PC)):
            print("Missing PC terms - fix required, NaNs found in PC")
            return
    
        print("No NaNs in PCs, continuing with reconstruction.")
        RC = np.zeros((N, M))
        for m in range(min(M, maxieg)):
            buf = np.outer(PC[:, m], RHO_maxieg[:, m].T)
            buf = np.flipud(buf)
            for n in range(N):
                RC[n, m] = np.mean(np.diagonal(buf, offset=-(N - M) + n))

    else:
         
        PC = np.dot(trajmat_ele_good, RHO)  # Principal components

        # Check for NaNs in the principal components
        if np.any(np.isnan(PC)):
            print("Missing PC terms - fix required, NaNs found in PC")
            return
        
        print("No NaNs in PCs, continuing with reconstruction.")

        # Calculate reconstructed components (RC)
        RC = np.zeros((N,M))
        for m in range(1,M+1):
            buf = np.outer(PC[:, m-1],RHO[:, m-1].T)  # Inverse projection
            buf = np.flipud(buf)  # Flip vertically

            for n in range(1,N+1):
                # Anti-diagonal averaging to reconstruct original components
                RC[n-1, m-1] = np.mean(np.diagonal(buf, offset=-(N - M ) + (n-1)))
    
    # Output structure containing results
    OutStruc = {
        'RC': RC,
        'LAMBDA': LAMBDA,
        'RHO': RHO,
        'num_nans': num_nans
    }

    print("Leaving local SSABasicCleanUCLAtoeplitz")
    return OutStruc

# %% Main code ---------------------------------------------------------------
# # Load data
directory = os.getcwd()  # Current directory    
input_dir ="To_SSA"  # Data directory

#print(f"Detected subfolders: {folders}")
file_names = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

for fldrs in file_names:
#fldrs = file_names[0]
    #fldrs = "CJTR_N"
    St_Names =fldrs[:10]
    #St_Names = St_Names.upper()
    Compnnt = fldrs[5]  # Component
    if Compnnt == '0':
        Comp = 'E'
    elif Compnnt == '1':
        Comp = 'N'
    elif Compnnt == '2':
        Comp = 'U'

    print(f"{St_Names} is processing" )

    filename = f"{St_Names}.dat"  # File match string
    filename_dir = os.path.join(input_dir, filename)  # Full file path

    comments = []
    data_lines = []

    with open(filename_dir, "r") as file:
        for line in file:
            if line.startswith("#"):
                comments.append(line.strip())  # Store comments without leading/trailing spaces
            else:
                data_lines.append(line)

    # Convert the numeric data to a NumPy array
    full_data = np.loadtxt(data_lines)


    n_data = full_data.shape[0]  # Number of data points
    data = full_data

    x_in = data[:, 3]  # Extract third column (mm) Observation - trend - jumps removed
    obser = data [:, 2] # Extract second column (mm) raw observation
    tyfrac = data[:, 1]  # Fractional year time
    tmjd = data[:, 0]  # Modified Julian Date

    dataindx = tmjd - tmjd[0] + 1  # Index of days with data

    # Normalize data
    x_std = np.std(x_in, ddof=1)
    x_mean = np.mean(x_in)
    x_in_mean_rem_norm_w_std = (x_in - x_mean) / x_std  # Mean removed and normalized

    # Create continuous time vector
    t = np.arange(tmjd[0], tmjd[-1] + 1) - tmjd[0] + 1
    n_cont_time = len(t)

    # Check for missing days
    if n_cont_time == n_data:
        print('No NaNs in input data')
    else:
        print(f'There are {n_cont_time - n_data} NaNs')



    # %% Create a copy of input vector with NaNs for missing values
    N = n_cont_time #number of continuous time steps

    x_in_with_NaNs = np.full((N,), np.nan)  # Equivalent to nan(n_cont_time,1) in MATLAB
    x_in_with_NaNs[dataindx.astype(int) - 1] = x_in  # Assign values at the correct indices

    # Create another NaN-filled array for normalized data
    x_in_mean_rem_norm_w_std_with_NaNs = np.full((N,), np.nan)
    x_in_mean_rem_norm_w_std_with_NaNs[dataindx.astype(int) - 1] = x_in_mean_rem_norm_w_std

    x = x_in_mean_rem_norm_w_std_with_NaNs  # This array has mean removed and normalized data

    #  Replace NaNs with zeros
    x_in_with_Zeros_mean_rem_norm_w_std = np.copy(x_in_mean_rem_norm_w_std_with_NaNs)
    x_in_with_Zeros_mean_rem_norm_w_std[np.isnan(x_in_with_Zeros_mean_rem_norm_w_std)] = 0

    X = x_in_with_Zeros_mean_rem_norm_w_std  # Final data with NaNs replaced by zeros


    # %% Plot gappy input time series (mean removed, normalized, NaNs for missing data)
    # Initialize figure number
    figno = 0
    #figname = f'{St_Names} Gappy input time series, Detrended, Jumps Removed'
    tFracYr_with_NaNs = np.full((N,), np.nan)  # Equivalent to nan(n_cont_time,1) in MATLAB
    tFracYr_with_NaNs[dataindx.astype(int) - 1] = tyfrac  # Assign values at the correct indices
    #plt.figure(figsize=(12, 9), dpi=300)
    #plt.plot(tFracYr_with_NaNs, x_in_with_NaNs, 'bo-', markersize=2, linewidth=0.3)  # 'bo-' -> blue circles with lines connecting
    #plt.xlabel('Time (days)')
    #plt.ylabel('Position')
    #plt.title(figname)
    Fig_output_dir = 'My_Figures'
    os.makedirs(Fig_output_dir, exist_ok=True)
    #fig_path = os.path.join(Fig_output_dir, f'1_{St_Names} Gappy input time series.png')
    #plt.savefig(f'{St_Names}Gappy_input_time_series_mean_removed_normalized_NaNs.png')
    #plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
    #plt.close()



    # %% Find gaps (NaNs)
    NaNindx = np.where(np.isnan(x))[0]  # Find indices of NaNs
    num_nans = len(NaNindx)  # Get number of NaNs
    frac_nans = num_nans / n_cont_time  # Fraction of NaNs

    # Check number of valid data points vs. total expected points
    n_data_ck = n_cont_time - num_nans

    def diff_ck(a1, a2, vari):
        """
        Check the difference between two arrays/scalars and print the max absolute difference.
        """
        del_val = np.max(np.abs(np.array(a1).flatten() - np.array(a2).flatten()))
        print(f"{vari}: {del_val}")

    diff_ck(n_data, n_data_ck, "Check: num data vs num gappy data in continuous")

    # %%  Embedding dimension
    # M is the number of columns in the trajectory matrix (embedding dimension)
    # Typically M is chosen as fix(n_data/2), but here it's manually set
    M = 910  # 2.5 years rounded to integer
    if len(tmjd) < 2 * M:
        print("Time series is shorter than 2*M, skipping...")
        continue
    # Number of rows in the trajectory matrix
    l = n_cont_time - M + 1  

    # %% Construct trajectory matrix
    colind = np.arange(1, M + 1)  # Column index vector (1:M in MATLAB)
    rowind = np.arange(0, l)[:, None]  # Row index vector (column vector)

    # Vectorized index matrix creation
    trajmatind = colind + rowind  # Equivalent to colind + rowind' in MATLAB

    # Initialize trajectory matrix with NaNs and fill values
    trajmat = np.full((l, M), np.nan)
    trajmat = x[trajmatind - 1]  # Adjust for Python's zero-based indexing

    # %% Toeplitz diagonal values calculation
    Cloop = np.zeros(M)  # Initialize Cloop array
    cntZeros_array = np.zeros(M)  # Store zero counts per J

    for J in range(M):
        Nterms = N - J  # Number of terms for this J
        cntZeros = 0  # Count zeros

        sum_val = 0  # Accumulator for sum
        for I in range(Nterms):
            newProdTerm = X[I] * X[I + J]
            if newProdTerm != 0:
                sum_val += newProdTerm
            else:
                cntZeros += 1

        ActualNterms = Nterms - cntZeros  # Adjust number of valid terms
        Cloop[J] = sum_val / ActualNterms if ActualNterms > 0 else 0  # Avoid division by zero

    # % Plot Toeplitz diagonal values
    #figname = f"{St_Names} Toeplitz diagonal values"
    #plt.figure(figsize=(12, 9), dpi=300)
    #plt.plot(Cloop, 'b-')
    #plt.xlabel('Lag')
    #plt.ylabel('Correlation')
    #plt.title(figname)
    #fig_path = os.path.join(Fig_output_dir, f'2_{St_Names}Toeplitz diagonal values.png')
    #plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
    #plt.close()


    # % Construct Toeplitz matrix
    Ctoep = toeplitz(Cloop)  # Generate Toeplitz matrix

    # % Plot Toeplitz lagged correlation matrix
    #figname = f"Toeplitz lagged correlation, gappy: Ntime={n_cont_time}, Ndata={n_data}, M={M}, L={l}"
    #plt.figure(figsize=(12, 9), dpi=300)
    #plt.imshow(Ctoep, aspect='auto', cmap='viridis')
    #plt.colorbar(label="Correlation")
    #plt.xlabel("Lag Index")
    #plt.ylabel("Lag Index")
    #plt.title(figname)
    #fig_path = os.path.join(Fig_output_dir, f'3_{St_Names} Toeplitz lagged correlation.png')
    #plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
    #plt.close()



    # %% Maximum eigenvalue
    maxeig=16

    # %% running the SSA code 
    OutStrucNaNs = SSABasicCleanUCLAtoeplitz(trajmat, toep=None, maxieg=maxeig)


    # %% FULL Reconstructed components (RC) from NaN-handled data
    RC = OutStrucNaNs['RC']

    # Insert NaNs into the RC matrix
    RCNaNs = RC.copy()
    if maxeig is None:
        RCNaNs[NaNindx, :] = np.nan  # Insert NaNs at specific indices

    # The above replaces the zeros (placeholders for NaNs) with NaNs again
        full_recons_NaNs = np.sum(RCNaNs, axis=1)  # Reconstruct full time series for NaN-handled data

    # Unnormalize the full reconstruction to compare with input (also process for velocity)
        x_in_recons_NaNs = full_recons_NaNs * x_std + x_mean

        # Plot full reconstruction
        figname = "gappy full reconstructions and data"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        plt.plot(tFracYr_with_NaNs, x_in_with_NaNs, label="input data")
        plt.plot(tFracYr_with_NaNs, x_in_recons_NaNs, label="full reconstruction")
        #plot_n(t, x, full_recons_NaNs, ['input data', 'full reconstruction'])
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(figname)
        fig_path = os.path.join(Fig_output_dir, f'4_{St_Names} Full_reconstruction.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'4_{St_Names} Full_reconstruction.png')

        figname = "Residual (Data - Reconstrction)"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        plt.plot(tFracYr_with_NaNs, x_in_with_NaNs - x_in_recons_NaNs, label="residuals")
        #plot_n(t, x - full_recons_NaNs, ['residuals'])
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(figname)
        fig_path = os.path.join(Fig_output_dir, f'5_{St_Names} Residuals.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'5_{St_Names} Residuals.png')


        maxeig_locl = 16

        
        # %% Plot eigenvalues
        figname = f"eigenvalues, max eig={maxeig_locl}"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        plt.semilogy(OutStrucNaNs['LAMBDA'], 'b+-')
        plt.semilogy(OutStrucNaNs['LAMBDA'][:maxeig_locl], 'ro-')
        plt.title(figname)
        fig_path = os.path.join(Fig_output_dir, f'6_{St_Names} Eigenvalues.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'6_{St_Names} Eigenvalues.png')

        figname = f"multiple part recon: Ntime={n_cont_time} Ndata={n_data}, M={M}, L={l}, num local eig={maxeig_locl}"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        # Call plot_n function
        plt.plot(tFracYr_with_NaNs, x, label="data")
        plt.plot(tFracYr_with_NaNs, np.sum(RCNaNs[:, 0:2], axis=1), label="Compnnt 1-2")
        plt.plot(tFracYr_with_NaNs, np.sum(RCNaNs[:, 2:4], axis=1), label="Compnnt 3-4")
        plt.plot(tFracYr_with_NaNs, np.sum(RCNaNs[:, 4:6], axis=1), label="Compnnt 5-6")
        plt.plot(tFracYr_with_NaNs, np.sum(RCNaNs[:, 6:maxeig_locl], axis=1), label=f'Compnnt 7-{maxeig_locl}')
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(f"{St_Names} SSA Reconstruction Components")  # Add title
        fig_path = os.path.join(Fig_output_dir, f'7_{St_Names} Multiple part reconstruction.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'7_{St_Names} Multiple part reconstruction.png')

        rc_maxEigenVal_lcl = np.sum(RCNaNs[:, :maxeig_locl], axis=1) * x_std + x_mean
        figname = "Reconstructed signals up to Local MAx Eigne Value"
        plt.figure(figname,figsize=(12, 9), dpi=300)
        figno = figno + 1
        # Call plot_n function
        plt.plot(tFracYr_with_NaNs, x_in_with_NaNs, label="data")
        plt.plot(tFracYr_with_NaNs, rc_maxEigenVal_lcl, label=f'Reconst to {maxeig_locl} Local Eigen value')
        #plot_n(t, x, rc_maxEigenVal_lcl, ['data', f'Reconst to {maxeig_locl} Local Eigen value'])
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(f"{St_Names} SSA Reconstruction Components")  # Add title
        fig_path = os.path.join(Fig_output_dir, f'8_{St_Names} Reconstructed signals Local Max Eigen Value.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'8_{St_Names} Reconstructed signals Local Max Eigen Value.png')

        RC_maxEigenVal = rc_maxEigenVal_lcl
    else:
        #NaNindx_maxeig = NaNindx[NaNindx < maxeig]  # Indices of NaNs within maxeig
        #RCNaNs[NaNindx_maxeig, :] = np.nan  # Insert NaNs at specific indices
        RCNaNs[NaNindx, :] = np.nan

            # %% Plot eigenvalues
        figname = f"{St_Names} eigenvalues, max eig={maxeig}"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        plt.semilogy(OutStrucNaNs['LAMBDA'], 'b+-')
        plt.semilogy(OutStrucNaNs['LAMBDA'][:maxeig], 'ro-')
        plt.title(figname)
        fig_path = os.path.join(Fig_output_dir, f'6_{St_Names} Eigenvalues.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'6_{St_Names} Eigenvalues.png')

        # %% Signals Reconstructed
        figname = f"{St_Names} multiple part recon: Ntime={n_cont_time} Ndata={n_data}, M={M}, L={l}, num eig={maxeig}"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        # Call plot_n function
        plt.plot(
            tFracYr_with_NaNs, x_in_with_NaNs, label='Data'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 0:maxeig], axis=1)* x_std + x_mean, label='Max Eigen Value Reconstruct'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 0:2], axis=1)* x_std + x_mean, label='Component 1-2'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 2:4], axis=1)* x_std + x_mean, label='Component 3-4'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 4:6], axis=1)* x_std + x_mean, label='Component 5-6'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 6:maxeig], axis=1)* x_std + x_mean, label=f'Component 7-{maxeig}'
        )
        plt.close()

        # %% Signals Reconstructed
        figname = f"{St_Names} 1st part recon: Ntime={n_cont_time}, num eig=1"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        # Call plot_n function
        plt.plot(
            tFracYr_with_NaNs, x_in_with_NaNs, label='Data'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 0:1], axis=1)* x_std + x_mean, label='Component 1'
        )
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(f"{St_Names} SSA Reconstruction 1st Components")  # Add title
        #plt.grid(True)  # Optional: add a grid for better visualization
        fig_path = os.path.join(Fig_output_dir, f'7_{St_Names} 1st part reconstruction.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        

                # %% Signals Reconstructed
        figname = f"{St_Names} 1st part recon: Ntime={n_cont_time}, num eig=2"
        figno = figno + 1
        plt.figure(figname,figsize=(12, 9), dpi=300)
        # Call plot_n function
        plt.plot(
            tFracYr_with_NaNs, x_in_with_NaNs, label='Data'
        )
        plt.plot(
            tFracYr_with_NaNs, np.sum(RCNaNs[:, 1:2], axis=1)* x_std + x_mean, label='Component 1'
        )
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(f"{St_Names} SSA Reconstruction 2nd Components")  # Add title
        #plt.grid(True)  # Optional: add a grid for better visualization
        fig_path = os.path.join(Fig_output_dir, f'7_{St_Names} 2nd part reconstruction.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'7_{St_Names} Multiple part reconstruction.png')

        # %% MAximum eigen value reconstruction
        # Compute rc by summing over the first maxeig columns
        rc_maxEigenVal = np.sum(RCNaNs[:, :maxeig], axis=1) * x_std + x_mean
        figname = f"{St_Names} Reconstructed signals up to MAx Eigne Value"
        plt.figure(figname,figsize=(12, 9), dpi=300)
        figno = figno + 1
        # Call plot_n function
        #plot_n(t, x, rc_maxEigenVal, ['data', f'Reconst to {maxeig} Eigen value'])
        plt.plot(
            tFracYr_with_NaNs, x_in_with_NaNs, label='Data'
        )
        plt.plot(
            tFracYr_with_NaNs, rc_maxEigenVal, '--',label=f'Reconst to {maxeig} Eigen value'
        )
        plt.legend()  # Add legend to distinguish the lines
        plt.xlabel("Time")  # Add X-axis label
        plt.ylabel("Amplitude")  # Add Y-axis label
        plt.title(f"{St_Names} Reconstructed signals Max Eigen Value: {maxeig}")  # Add title
        fig_path = os.path.join(Fig_output_dir, f'8_{St_Names} Reconstructed signals Max Eigen Value {maxeig}.png')
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
        plt.close()

        #plt.savefig(subfolder, f'8_{St_Names} Reconstructed signals Max Eigen Value {maxeig}.png')

        RC_maxEigenVal = rc_maxEigenVal  # Reconstructed components up to maxeig
    # %%Lomb Scargle 
    # Compute FFT

    '''
    What Does Lomb-Scargle Represent?
    The Lomb-Scargle Periodogram is a method for estimating the power spectrum of a signal, especially when data points are irregularly spaced in time. It is commonly used in astronomy, geophysics, and signal processing to detect periodic signals in unevenly sampled data.

    Key Concepts
    Power Spectrum Estimation

    It calculates how much power (variance) in the data is associated with different frequencies.
    Peaks in the periodogram indicate dominant frequencies (periodic components) in the data.
    '''

    residual_RC = x - RC_maxEigenVal

    # Assume t and rc_maxEigenVal are your original time and signal arrays
    rc_mask = ~np.isnan(RC_maxEigenVal)  # Create a mask for valid (non-NaN) values
    residual_mask = ~np.isnan(residual_RC)  # Create a mask for valid (non-NaN) values

    # Apply mask to remove NaN values from both arrays
    t_rc_clean = t[rc_mask]
    t_res_clean = t[residual_mask]

    rc_clean = RC_maxEigenVal[rc_mask]
    res_clean = residual_RC[residual_mask]

    # from astropy.timeseries import LombScargle
    # frequency, power = LombScargle(t_clean, rc_clean).autopower(minimum_frequency=0.1, maximum_frequency=0.5,samples_per_peak=10)


    #from scipy.signal import lombscargle
    from astropy.timeseries import LombScargle
    w = np.linspace(0.01, 0.5, 250)
    min_freq = 0.1
    frequency_RC, pgram_powerRC = LombScargle(t_rc_clean, rc_clean, normalization='psd').autopower(minimum_frequency=min_freq, maximum_frequency=100)
    frequency_Noise, pgram_powerNoise = LombScargle(t_res_clean, res_clean, normalization='psd').autopower(minimum_frequency=min_freq, maximum_frequency=100)

    #pgram_power = lombscargle(t_rc_clean, rc_clean, w, normalize=True)

    #Res_RC_power = lombscargle(t_res_clean, res_clean, w, normalize=True)
    #power_ratio = np.sum(pgram_power) / (np.sum(pgram_power) + np.sum(Res_RC_power))
    #If power_ratio is close to 1, the selected components capture most cyclic behavior.
    #If power_ratio is small, cyclic components are still present in the residual signal.

    #figname = "Lomb Scargle for the reconstructed Signal"
    #figno = figno + 1
    plt.figure(figsize=(12, 9), dpi=300)
    plt.loglog(frequency_RC, pgram_powerRC, 'r',markersize=1, label='cyclics')
    plt.loglog(frequency_Noise, pgram_powerNoise, 'b.-',markersize=1,label='noise')
    #plt.title(f"Lomb-Scargle Power Spectrum for {maxeig}\nPower Ratio: {power_ratio:.2f}")
    plt.xlabel('freq')
    plt.ylabel('Amp')
    plt.title(f'{St_Names} Lomb Scargle Amplitude Spectrum')
    plt.grid()
    plt.legend()
    fig_path = os.path.join(Fig_output_dir, f'9_{St_Names} Lomb Scargle Amplitude Spectrum.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')  # Save the figure with tight layout
    plt.close()
    #plt.savefig(subfolder, f'9_{St_Names} Lomb Scargle Amplitude Spectrum.png')

    # %% Saving the data into a file


    valid_indices = ~np.isnan(RC_maxEigenVal)  # Mask for valid (non-NaN) values
    RC_maxEigenVal_filtered = RC_maxEigenVal[valid_indices]

    # Define the output file name and path
    output_filename = f"{St_Names}_Cyclics.mom"
    SSA_output_dir = 'SSA_Cyclics'
    os.makedirs(SSA_output_dir, exist_ok=True)  # ✅ Only the directory
    output_path = os.path.join(SSA_output_dir, output_filename)

    # Open a file and write the data
    with open(output_path, "w", newline="") as file:
        for t, r in zip(tmjd, RC_maxEigenVal_filtered):
            file.write(f"{t} {r:.6f}\n")
    # Save the data into a file
            
            
    sampling_period = "# sampling period 1.0"
    output_filename_mom = f"{St_Names}_{Comp}.mom"
    output_dir_mom = os.path.join(directory, "..", "4th_Velocity_HectorP_SSACycleRmved")
    subfolder_mom = os.path.join(output_dir_mom,"obs_files")
    os.makedirs(subfolder_mom, exist_ok=True)  # Create the subfolder if it doesn't exist 
    with open(os.path.join(subfolder_mom,output_filename_mom), "w", newline="") as file_mom:
            # Write the sampling period comment
            file_mom.write(sampling_period + "\n")
            # Write the stored comment lines
            if comments and comments[0].startswith("# exp"):
                comments[0] = comments[0].replace("# exp", "# log", 1)
            for comment in comments:
                file_mom.write(comment + "\n")
            for ti, obs in zip(tmjd, obser-RC_maxEigenVal_filtered):
                file_mom.write(f"{ti} {obs:.6f}\n")  # Format r to 6 decimal places
