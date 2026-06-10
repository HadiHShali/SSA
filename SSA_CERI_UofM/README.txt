SSA_Schoellhamer_UofM_V3.py
============================
Singular Spectrum Analysis (SSA) for GNSS time series.
Reads daily .dat files, decomposes each series into reconstructed components (RCs),
saves SSA cyclic outputs and Hector-format residual .mom files, and optionally
generates diagnostic figures.


REQUIREMENTS
------------
    numpy  scipy  matplotlib  astropy  pandas


USAGE
-----
    python SSA_Schoellhamer_UofM_V3.py [options]

Minimal:
    python SSA_Schoellhamer_UofM_V3.py --MaxEig 16 --plot-reconstruct --plot-EigValue --plot-lomb

Full:
    python SSA_Schoellhamer_UofM_V3.py
        --input  ToSSA_Input/input_TS/
        --output-ts  SSA_OutPut/output_TS/
        --output-fig SSA_OutPut/output_Figures/
        --embedding 910 --MaxEig 16
        --components "1-2,3-6,7-16"
        --plot-all


ARGUMENTS
---------
I/O:
  --input          Directory of .dat input files          (default: ToSSA_Input/input_TS/)
  --output-ts      SSA cyclic time-series outputs         (default: SSA_OutPut/output_TS/)
  --output-fig     Figure root, one sub-folder/station    (default: SSA_OutPut/output_Figures/)

Processing:
  --embedding/-M   Embedding dimension in days (~2.5yr)   (default: 910)
  --MaxEig         Eigenvalues used in reconstruction     (default: 16)
  --components     RC groupings for multi-panel plot      (default: "1-2,3-4,5-6")
  --toeplitz       Use Toeplitz covariance matrix         (flag, off by default)
  --stations       Process only listed station names      (default: all)
  --dry-run        Print stations to process then exit    (flag)

Lomb-Scargle PSD:
  --min-freq       Lower frequency bound (cycles/year)    (default: 0.3)
  --max-freq       Upper frequency bound (cycles/year)    (default: 100.0)

Plot flags:
  --plot-all         Enable every plot
  --plot-gappy       Gappy input time series
  --plot-toeplitz    Toeplitz diagonal and correlation matrix
  --plot-EigValue    Eigenvalue spectrum (all=blue, selected=red)
  --plot-reconstruct Summed reconstruction vs data + RC group panel
  --plot-individuals One figure per RC component (1 to MaxEig)
  --plot-residuals   Residuals (data minus reconstruction)
  --plot-lomb        Lomb-Scargle PSD of reconstruction and residuals


--components SYNTAX
-------------------
Comma-separated 1-based ranges, e.g.:

    --components "1-2,3-6,7-16"

    1-2   -> sums RC 1 and 2,        label: RC 1-2
    3-6   -> sums RC 3 through 6,    label: RC 3-6
    7-16  -> sums RC 7 through 16,   label: RC 7-16

If a group's end exceeds --MaxEig it is clipped to MaxEig.
The label reflects the actual range used.


OUTPUT STRUCTURE
----------------
    SSA_OutPut/
    +-- output_TS/
    |   +-- STATIONNAME_Cyclics.mom         (SSA cyclic time series)
    +-- output_Figures/
        +-- STATIONNAME/
            +-- 01_gappy_input.png
            +-- 02_toeplitz.png
            +-- 03_eigenvalues.png
            +-- 04_reconstruction.png
            +-- 05_rc_groups.png
            +-- 06_rc_01.png ... 06_rc_16.png
            +-- 07_residuals.png
            +-- 08_lomb_scargle.png

    ../4th_Velocity_HectorP_SSACycleRmved/obs_files/
        +-- STATIONNAME_E.mom / _N.mom / _U.mom   (Hector-format residuals)


INPUT FILE FORMAT
-----------------
Plain-text .dat files; comment lines start with #; data columns:

    MJD   FractionalYear   Observation_mm   Detrended_mm

File name pattern: STATIONNAME.dat  (10-character station name).
Character at index 5 encodes the component: 0=E, 1=N, 2=U.
