# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:37:44 2024

@author: GeodesyLab
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz

# Clear variables and close all figures
plt.close('all')

# Set the precision for displaying numbers (similar to MATLAB's format shortG)
np.set_printoptions(precision=7, suppress=True)

# Initialize variables
figno = 0

# Process Hector file from Hadi
filename = '1LSU_1_SSA.dat'  # Input data file with gaps (MJD, fracyr, posn; time not continuous)


## --------------------------------------------------------------------------##
## -------------------------------Functions----------------------------------##
## --------------------------------------------------------------------------##
# (1)
# Function to read the data (assuming a simple whitespace-delimited text file)
def read_matrix(file):
    try:
        return np.loadtxt(file)
    except Exception as e:
        print(f"Error reading the file: {e}")
        return None
    
    
    
    
    
    
    
    
    
    
    
    
## --------------------------------------------------------------------------##
## -------------------------------Input Data---------------------------------##
## --------------------------------------------------------------------------##
    
# Read the data
full_data = read_matrix(filename)

# Note: There's a big natural gap in data from day 4847 (last day) to 5085 (first day back), missing 238 days in between.

# This is the number of data points (not the same as the number of days from start to finish due to gaps)
n_data = full_data.shape[0]

#data = full_data : Changes to data will also affect full_data.
data = full_data.copy()

# Extract Columns
x_in = data[:,2] # Ordered but not uniformly spaced, column vector
tyfrac = data[:,1] # Time in fraction of year (good for plotting, not uniformly spaced)
tmjd = data[:,0]  # Time of data in modified Julian days, good for indexing (not uniformly spaced)

# Compute indices of days with data
dataindx = tmjd - tmjd[0] + 1  # Indices of days with data (will line up with continuous t later)

# Remove mean and normalize
x_mean = np.mean(x_in)
x_std = np.std(x_in)
x_in_mean_rem_norm_w_std = (x_in - x_mean) / x_std  # Data for SSA, mean removed and normalized by SD

# Make a continuous time vector that starts from 1
t = np.arange(tmjd[0], tmjd[-1] + 1) - tmjd[0] + 1  # Continuous time vector
n_cont_time = len(t)  # Length of continuous time vector

# Check if there are gaps
if n_cont_time == len(x_in):
    print("No NaNs in input data")
else:
    print(f"There are {n_cont_time - len(x_in)} NaNs")
    
N = n_cont_time  # Continuous time vector length

# Create a copy of the input vector with NaNs for continuous time length
x_in_with_NaNs = np.full((n_cont_time,), np.nan)
x_in_with_NaNs[dataindx.astype(int) - 1] = x_in  # Place data values where indices match

# Create a normalized version with NaNs
x_in_mean_rem_norm_w_std_with_NaNs = np.full((n_cont_time,), np.nan)
x_in_mean_rem_norm_w_std_with_NaNs[dataindx.astype(int) - 1] = x_in_mean_rem_norm_w_std

# Assign to `x`, this has mean removed and normalized to standard deviation
x = x_in_mean_rem_norm_w_std_with_NaNs

# Replace NaNs with zeros in the normalized series
x_in_with_Zeros_mean_rem_norm_w_std = np.copy(x_in_mean_rem_norm_w_std_with_NaNs)
x_in_with_Zeros_mean_rem_norm_w_std[np.isnan(x_in_with_Zeros_mean_rem_norm_w_std)] = 0
X = x_in_with_Zeros_mean_rem_norm_w_std  # Final version with NaNs replaced by zeros

## --------------------------Plot the Input Data-----------------------------##
figname = 'Gappy input time series, mean removed and normalized, NaNs for missing data'

# Create a new figure
plt.figure()
plt.title(figname)

# Plot the data
plt.plot(t, x, 'bo-', label='Normalized Data', markersize=2, linewidth=0.8)

# Add labels to the axes
plt.xlabel('Time (days)')
plt.ylabel('Position')

# Display the legend
plt.legend()

# Show the figure
plt.show()
## --------------------------------------------------------------------------##


## -----------Check the NaNs --------------------##
# Find indices of NaNs
NaNindx = np.where(np.isnan(x))[0]  # Indices of NaNs
num_nans = len(NaNindx)  # Number of NaNs

# Fraction of NaNs in the continuous time series
frac_nans = num_nans / n_cont_time  # Fraction of NaNs

# Number of non-NaN data points
n_data_ck = n_cont_time - num_nans

# Check if the number of original data points matches the count after accounting for NaNs
assert n_data == n_data_ck, f"Check failed: n_data ({n_data}) != n_data_ck ({n_data_ck})"


## --------------------------------------------------------------------------##
## ----------------Trajectory Matrix / Correlation Matrix--------------------##
## --------------------------------------------------------------------------##

# Set the embedding dimension (M), which determines the number of columns in the trajectory matrix
M = 1278  # Embedding dimension set to 1278, representing approximately 3.5 years of data points

# Calculate the number of rows (l) in the trajectory matrix
l = n_cont_time - M + 1  # Number of rows depends on the total data length and embedding dimension

# Create the column index vector (1 to M) for the trajectory matrix
colind = np.arange(0, M)  # Column indices range from 1 to M

# Create the row index vector (0 to l-1) for the trajectory matrix
rowind = np.arange(0, l)  # Row indices range from 0 to l-1

# Generate trajectory matrix indices by combining column and row indices
trajmatind = colind + rowind[:, None]  # This vectorized operation combines column and row indices

# Initialize the trajectory matrix with NaN values
trajmat = np.full((l, M), np.nan)  # Create a matrix of size (l, M) filled with NaNs

# Fill the trajectory matrix with data from the normalized time series (x)
trajmat = x[trajmatind-1]  # Populate the matrix using the values from x based on the calculated indices

## ----------------------Plot the Trajectory Matrix--------------------------##
plt.figure(figsize=(10, 6))
plt.imshow(trajmat, aspect='auto', cmap='viridis', interpolation='nearest')
plt.colorbar(label='Trajectory Values')  # Add color bar to show the value scale
plt.title('Trajectory Matrix')
plt.xlabel('Embedding Dimension')
plt.ylabel('Time Steps')
plt.show()
## --------------------------------------------------------------------------##

#do with toeplitz matrix

Cloop = np.zeros(M)  # Initialize Cloop array

# Loop to compute the diagonal values of the Toeplitz matrix
for J in range(M):
    Jindx = J  # Index in Python starts at 0, so no need to add 1
    Cloop[Jindx] = 0
    Nterms = N - J
    cntZeros = 0
    for I in range(Nterms):
        newProdTerm = X[I] * X[I + J]
        if newProdTerm != 0:
            Cloop[Jindx] += newProdTerm
        else:
            cntZeros += 1
    ActualNterms = Nterms - cntZeros
    Cloop[Jindx] = Cloop[Jindx] / ActualNterms

## ----------------------Plot the Cloop--------------------------##
figname = f"Toeplitz Diagonal Values (Gappy): Ntime={n_cont_time}, Ndata={n_data}, M={M}, L={l}"
plt.figure(figsize=(10, 6))
plt.plot(Cloop)
plt.title(figname)
plt.xlabel('Lag (J)')
plt.ylabel('Correlation Value')
plt.show()
## --------------------------------------------------------------------------##

# Create the Toeplitz matrix
Ctoep = toeplitz(Cloop)

## ----------------------Plot the Toeplitz Matrix--------------------------##
figname = f"Toeplitz Lagged Correlation (Gappy): Ntime={n_cont_time}, Ndata={n_data}, M={M}, L={l}"
plt.figure(figsize=(10, 6))
plt.imshow(Ctoep, aspect='auto', cmap='viridis')
plt.colorbar(label='Correlation Value')
plt.title(figname)
plt.xlabel('Lag (J)')
plt.ylabel('Lag (J)')
plt.show()
## --------------------------------------------------------------------------##


## --------------------------------------------------------------------------##
## ------------------------------SSA Analysis--------------------------------##
## --------------------------------------------------------------------------##
#OutStrucNaNs = SSABasicCleanUCLAtoeplitz(trajmat)

























