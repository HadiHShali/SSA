# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:21:31 2025

@author: GeodesyLab
"""
# %matplotlib qt

# %% Clear workspace and setup
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz  # For Toeplitz matrix generation

# %% Functions
#plot_n function (Function to plot multiple datasets with legends and optional residuals)
def plot_n(*args):
    t = args[0]
    legendlabel = args[-1]
    numplots = len(legendlabel)
    n_vec_in = len(args) - 2

    if numplots == 2:
        symcolr = ['c', 'r--']
    else:
        symcolr = ['c', 'r', 'b', 'k', 'm']

    #plt.hold(True)
    y = np.zeros((len(t), n_vec_in))
    for ydata_n in range(n_vec_in):
        y[:, ydata_n] = args[ydata_n + 1]
        plt.plot(t, y[:, ydata_n], symcolr[ydata_n])
    
    plt.xlabel('time (days)')
    plt.ylabel('position')

    numres = numplots - n_vec_in
    scl = 10
    offset = -3
    for resno in range(numres):
        res = y[:, 0] - y[:, resno + 1]
        maxres = np.max(np.abs(res))
        plt.plot(t, scl * (y[:, 0] - y[:, resno + 1]) + offset, symcolr[resno + 2])

    plt.legend(legendlabel, loc='lower right')


# SSA function
# % SSA Base Function
def SSABasicCleanUCLAtoeplitz(trajmat, toep=None):
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
    if num_nans > 0:
        print("NaNs detected in principal components, need to handle them.")
    
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
# %% Load data
filename = '1NSU_0_SSA.dat'  # Input data file

full_data = np.loadtxt(filename)

# %% Data processing
n_data = full_data.shape[0]  # Number of data points
data = full_data

x_in = data[:, 2]  # Extract third column
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
figname = 'Gappy input time series, mean removed and normalized, NaNs for missing data'

plt.figure()
plt.plot(t, x, 'bo-', markersize=2, linewidth=0.3)  # 'bo-' -> blue circles with lines connecting
plt.xlabel('Time (days)')
plt.ylabel('Position')
plt.title(figname)
plt.show()

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
M = 1278  # 3.5 years rounded to integer

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
figname = f"Toeplitz diagonal values double loop - mine 1, gappy: Ntime={n_cont_time}, Ndata={n_data}, M={M}, L={l}"
plt.figure()
plt.plot(Cloop, 'b-')
plt.xlabel('Lag')
plt.ylabel('Correlation')
plt.title(figname)
plt.show()

# % Construct Toeplitz matrix
Ctoep = toeplitz(Cloop)  # Generate Toeplitz matrix

# % Plot Toeplitz lagged correlation matrix
figname = f"Toeplitz lagged correlation, gappy: Ntime={n_cont_time}, Ndata={n_data}, M={M}, L={l}"
plt.figure()
plt.imshow(Ctoep, aspect='auto', cmap='viridis')
plt.colorbar(label="Correlation")
plt.xlabel("Lag Index")
plt.ylabel("Lag Index")
plt.title(figname)
plt.show()

# %% running the SSA code 
OutStrucNaNs = SSABasicCleanUCLAtoeplitz(trajmat, toep=None)

maxeig=50;

# %% Plot eigenvalues
figname = f"eigenvalues, max eig={maxeig}"
figno = figno + 1
fig_handle = plt.figure(figname)
plt.semilogy(OutStrucNaNs['LAMBDA'], 'b+-')
plt.semilogy(OutStrucNaNs['LAMBDA'][:maxeig], 'ro-')
plt.title(figname)
plt.show()


# %% FULL Reconstructed components (RC) from NaN-handled data
RC = OutStrucNaNs['RC']

# Insert NaNs into the RC matrix
RCNaNs = RC.copy()
RCNaNs[NaNindx, :] = np.nan  # Insert NaNs at specific indices

# The above replaces the zeros (placeholders for NaNs) with NaNs again
full_recons_NaNs = np.sum(RCNaNs, axis=1)  # Reconstruct full time series for NaN-handled data

# Unnormalize the full reconstruction to compare with input (also process for velocity)
x_in_recons_NaNs = full_recons_NaNs * x_std + x_mean

# Plot full reconstruction
figname = "gappy full reconstructions and data"
figno = figno + 1
fig_handle = plt.figure(figname)
plot_n(t, x, full_recons_NaNs, ['input data', 'full reconstruction'])
plt.title(figname)
plt.show()

figname = "Residual (Data - Reconstrction)"
figno = figno + 1
fig_handle = plt.figure(figname)
plot_n(t, x - full_recons_NaNs, ['residuals'])
plt.title(figname)
plt.show()

# %% Signals Reconstructed
figname = f"multiple part recon: Ntime={n_cont_time} Ndata={n_data}, M={M}, L={l}, num eig={maxeig}"
figno = figno + 1
fig_handle = plt.figure(figname)
# Call plot_n function
plot_n(
    t, x, 
    np.sum(RCNaNs[:, 0:2], axis=1),
    np.sum(RCNaNs[:, 2:4], axis=1),
    np.sum(RCNaNs[:, 4:6], axis=1),
    np.sum(RCNaNs[:, 6:maxeig], axis=1),
    ['data', 'Comp 1-2', 'Comp 3-4', 'Comp 5-6', f'Comp 7-{maxeig}']
)

plt.show()

# %% MAximum eigen value reconstruction
# Compute rc by summing over the first maxeig columns
rc_maxEigenVal = np.sum(RCNaNs[:, :maxeig], axis=1)
figname = "Reconstructed signals up to MAx Eigne Value"
fig_handle = plt.figure(figname)
figno = figno + 1
# Call plot_n function
plot_n(t, x, rc_maxEigenVal, ['data', f'Reconst to {maxeig} Eigen value'])
plt.show()


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

residual_RC = x - rc_maxEigenVal

# Assume t and rc_maxEigenVal are your original time and signal arrays
rc_mask = ~np.isnan(rc_maxEigenVal)  # Create a mask for valid (non-NaN) values
residual_mask = ~np.isnan(residual_RC)  # Create a mask for valid (non-NaN) values

# Apply mask to remove NaN values from both arrays
t_rc_clean = t[rc_mask]
t_res_clean = t[residual_mask]

rc_clean = rc_maxEigenVal[rc_mask]
res_clean = residual_RC[residual_mask]

# from astropy.timeseries import LombScargle
# frequency, power = LombScargle(t_clean, rc_clean).autopower(minimum_frequency=0.1, maximum_frequency=0.5,samples_per_peak=10)


from scipy.signal import lombscargle
w = np.linspace(0.01, 0.5, 250)
pgram_power = lombscargle(t_rc_clean, rc_clean, w, normalize=True)

Res_RC_power = lombscargle(t_res_clean, res_clean, w, normalize=True)
#power_ratio = np.sum(pgram_power) / (np.sum(pgram_power) + np.sum(Res_RC_power))
#If power_ratio is close to 1, the selected components capture most cyclic behavior.
#If power_ratio is small, cyclic components are still present in the residual signal.

figname = "Lomb Scargle for the reconstructed Signal"
figno = figno + 1
fig_handle = plt.figure(figname)
plt.loglog(w, pgram_power, 'r',markersize=1, label='cyclics')
plt.loglog(w, Res_RC_power, 'b.-',markersize=1,label='noise')
#plt.title(f"Lomb-Scargle Power Spectrum for {maxeig}\nPower Ratio: {power_ratio:.2f}")
plt.xlabel('freq')
plt.ylabel('Amp')
plt.grid()
plt.legend()
plt.show()

