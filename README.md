# Singular Spectrum Analysis (SSA) for GNSS Time Series

> **A Python implementation of Singular Spectrum Analysis for decomposing GNSS time series into their underlying frequency components, with native support for missing data (gaps).**

---

## What This Project Does

Singular Spectrum Analysis (SSA) is a non-parametric spectral technique that **decomposes a time series into its constituent oscillatory components** — without assuming any specific functional form for the signals.

In simple terms:
-  Take a noisy GNSS time series (it's better to first remove the jump(s) and trend from the timeseries)
-  Decompose it into a sum of frequencies (seasonal cycles, oscillations, noise, etc.)
-  Reconstruct using all components → recover the **exact original signal** (residual ≈ 0)
-  Or reconstruct using only selected frequencies → isolate phenomena of interest

This makes SSA especially powerful for analyzing **GNSS station displacement records** where one needs to separate tectonic signals from hydrological loading, seasonal motion, and noise.

### Why Gaps Matter

Real GNSS time series almost always contain gaps due to equipment failures, weather, or maintenance. Classical SSA cannot handle gaps directly. This implementation adapts the **Schoellhamer (2001)** gap-tolerant SSA algorithm — making it usable on real-world GNSS data without interpolation.

![GPS Station vs River Gauge Comparison](https://github.com/user-attachments/assets/443cc1bf-cbcf-4fe8-99dd-1004966ed809)

*Example: SSA decomposition of GPS station LCHS displacement compared against nearby river gauge MS116 — reveals common hydrologically-driven signals.*

---

## Workflow

```
1- Raw GNSS Time Series (with gaps)
2- Hectorp or any other packages that you know (to remove trend and jumps)
3- SSA (handle gaps, decompose)
4- Frequency Components + Reconstructions

```

**Pre-processing** — Remove offsets, trends, and jumps using the [Hectorp](https://gitlab.com/machielsimonbos/hectorp) package. Or any other pacakges you know. **We only need an input time series without Jums and Trend.** 
**Hectorp Package** - I show you how to detrend and remove the jumps using hectorp module in python. 
**Gap-tolerant SSA** — Apply our modified Schoellhamer (2001) algorithm to decompose the cleaned signal
**Component analysis** — Identify and group reconstructed components (RCs) by frequency
**Signal reconstruction** — Recombine selected components to isolate signals of interest

---

## Repository Structure

```
SSA/
│
├── 1_Papers/                      				# Reference papers (the theory)
│   └── Schoellhamer_2001.pdf      				#   - Base method for gap-tolerant SSA
│
├── 2_UCLA_Code/                   				# Reference MATLAB implementation
│   └── ...                        				#   - Standard SSA tutorial (no gap support)
│
├── 3_Schoellhamer_Code/           				# Adapted MATLAB implementation
│   ├── original/                  				#   - Schoellhamer's original code
│   └── modified/                  				#   - Our modifications for GNSS workflows
│
├── 4_Main_SSA_CERI/               				# Main Python pipeline 
│   ├── SSA_Schoellhamer_UofM_V3.py             #   - Core SSA implementation
│   ├── ToSSA_Input/                   				#   - Worked example with sample data
│   │   ├── input_TS/
│   │   ├── NevGeoLab/
│   ├── SSA_OutPut/
│   │   └── output_Figures/
│   │   └── output_TS/
│   └── README.txt                  				#   - Detailed usage guide
│
└── README.md                      				# This file
```

---

## Quick Start

### Prerequisites
I'm using windows 11 CMD here. 
 1. Verify Python version (3.9+)
        python3 --version

 2. Download the repository and make the required steps:
 2.1. Clone the repo:
        git clone https://github.com/HadiHShali/SSA.git
 2.2. Go to the folder:
        cd SSA

 3. Install required packages inside the environment
        python.exe -m pip install --upgrade pip
        pip install numpy pandas requests scipy matplotlib hectorp

### Run the Example (in Windows CMD)
### 1- your own data without trend and jumps
Read the README.txt file inside the SSA_CERI_UofM folder. 
Put your data inside the ToSSA_Input\input_TS 
Run the code as instructed in the REAME.txt file. 

### 2- Nevada Geodetic Labrotary data (use Hector package to detrend and remove the jumps)
1- Download the data from Nevada Geodetic Lab (tenv3 files):
(the following peice of code give you all the stations in North American plate loctaed in 35<lat<36 and -90<lon<-89 and have minimum 5 years of data)

python download_gnss_stations.py --ref NA --min-lat 35.0 --max-lat 36.0 --min-lon -90 --max-lon -89 --min-duration 5

It dowloads the .tenv3 files and saves them in tenv3 folder

2- Read the tenv3 files and make directories for each components withing the Stns_Dir

python tenv3_reader_compDir.py

3- While the virtual environmet is activated:

### 1- your own raw data (use Hector package to detrend and remove the jumps)
