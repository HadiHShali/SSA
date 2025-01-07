# Singular Spectrum Analysis (SSA)
This is the repository for SSA which include different sections as follow:

# 1- Papers
This folder includes all the papers that we used to understand, implement, and test the SSA method for the data containing gaps.
 
The main and base paper is "Singular spectrum analysis for time series with missing data" by David H. Schoellhamer, U.S. Geological Survey, Sacramento, California. 
GEOPHYSICAL RESEARCH LETTERS, VOL. 0, NO. 0, PAGES 0-0, M 0, 2001.

# 2- UCLA code
SSA and M-SSA tutorial in Matlab, very suitable for understanding the steps of Singular Spectrum Analysis. 

Links:

https://www.mathworks.com/matlabcentral/fileexchange/58967-singular-spectrum-analysis-beginners-guide

https://www.mathworks.com/matlabcentral/fileexchange/58968-multichannel-singular-spectrum-analysis-beginners-guide

This Matlab tutorial demonstrates step by step the univariate as well as multivariate singular spectrum analysis. The steps are almost similar to those of a singular spectrum analysis.

References

[1] Groth, A., and M. Ghil, 2015: Monte Carlo Singular Spectrum Analysis (SSA) revisited: Detecting oscillator clusters in multivariate datasets, Journal of Climate, 28, 7873-7893.

Find more research on Singular Spectrum Analysis and related topics at https://dept.atmos.ucla.edu/tcd

# 3- Schoellhamer code
The original Schoellhamer code as well as its modified version to use for our data are included in this section. 

The original code consists of two functions, one for SSA with Nan and the other for reconstructed components. 

In our modified version, we just made some modifications as well as a function that reads the output data from Hectorp pacakge. In details, we have removed the jumps and trend using the hectorp package and feed into the code (Schoellhamer modified version) to do the SSA. 

![Scholhamer test](https://github.com/user-attachments/assets/8ea3128a-e0b8-4285-909b-cb6fa0ae07b2)

Reference(s):

[1] Schoellhamer, David H. "Singular spectrum analysis for time series with missing data." Geophysical research letters 28.16 (2001): 3187-3190.

# 4- Main SSA of CERI, UofM
