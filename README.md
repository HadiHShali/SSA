<div align="center">

# Gap-Tolerant Singular Spectrum Analysis for GNSS Time Series

**A Python toolkit for decomposing GNSS displacement records into their underlying frequency components — with first-class support for missing data (gaps), no interpolation required.**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)]()
[![Method](https://img.shields.io/badge/Method-Schoellhamer%202001-success.svg)]()
[![Status](https://img.shields.io/badge/Status-Active-brightgreen.svg)]()
<!-- Add a license badge once you've chosen one, e.g. MIT: -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) -->

</div>

---

## Overview

Singular Spectrum Analysis (SSA) is a **non-parametric** spectral technique that decomposes a time series into a set of interpretable, data-driven components — trends, quasi-periodic cycles, and noise — **without assuming any functional form** for the underlying signals.

This makes it a natural fit for GNSS station displacement records, where the goal is often to **separate tectonic signal from hydrological loading, seasonal motion, and colored noise**. Because it learns structure directly from the data, SSA is well suited to signals that parametric models struggle to capture cleanly.

The key differentiator of this implementation is **robustness to real-world data**: classical SSA breaks down in the presence of gaps, but GNSS records almost always contain them. This toolkit adapts the gap-tolerant algorithm of **Schoellhamer (2001)** so it works directly on operational GNSS data.

<div align="center">

![GPS Station vs River Gauge Comparison](https://github.com/user-attachments/assets/443cc1bf-cbcf-4fe8-99dd-1004966ed809)

*SSA decomposition of GPS station **LCHS** displacement compared against nearby river gauge **MS116** — isolating a common, hydrologically-driven signal from the raw record.*

</div>

---

## Key Features

- **Gap-tolerant by design** — handles missing epochs natively using a modified Schoellhamer (2001) lagged-covariance approach. No gap-filling or interpolation needed.
- **Exact reconstruction** — recombine all components to recover the original signal to numerical precision (residual ≈ 0), or select a subset to isolate the phenomena you care about.
- **Non-parametric decomposition** — no need to pre-specify periods, harmonics, or a functional model; the method extracts structure directly from the data.
- **Built for GNSS workflows** — integrates cleanly with [Hectorp](https://gitlab.com/machielsimonbos/hectorp) for offset/trend removal and with Nevada Geodetic Laboratory data products.
- **Reproducible, end-to-end pipeline** — from raw `.tenv3` download through detrending, decomposition, and reconstruction, with a worked example included.

---

## Why Gaps Matter

Real GNSS time series are rarely continuous. Equipment failures, power and communication outages, weather, and routine maintenance all introduce gaps. Classical SSA cannot ingest these series without interpolating across the missing intervals — which injects artifacts and biases the recovered components.

This implementation instead estimates the lagged-covariance structure from the **available data only**, following Schoellhamer (2001). The result is a decomposition that reflects the true signal rather than the interpolation scheme — a prerequisite for trustworthy analysis on operational data.

---

## How It Works

```
  Raw GNSS Time Series (with gaps)
              │
              ▼
  Pre-processing  ──  remove offsets, jumps, and trend (Hectorp or equivalent)
              │
              ▼
  Gap-Tolerant SSA  ──  decompose into data-driven components
              │
              ▼
  Component Analysis  ──  group reconstructed components (RCs) by frequency
              │
              ▼
  Signal Reconstruction  ──  recombine selected components to isolate signals
```

**1 · Pre-processing** — Remove offsets, jumps, and the linear trend so the input to SSA is a clean, stationary-ish residual. We use [Hectorp](https://gitlab.com/machielsimonbos/hectorp) and show how below, but any equivalent tool works. *SSA only needs an input series with jumps and trend already removed.*

**2 · Gap-tolerant SSA** — Apply the modified Schoellhamer (2001) algorithm to decompose the cleaned signal.

**3 · Component analysis** — Identify and group the reconstructed components (RCs) by frequency.

**4 · Signal reconstruction** — Recombine the selected components to isolate the signal of interest (e.g. seasonal cycles, hydrological loading, or residual tectonic motion).

---

## Repository Structure

```
SSA/
│
├── 1_Papers/                          # Reference papers (the theory)
│   └── Schoellhamer_2001.pdf          #   Base method for gap-tolerant SSA
│
├── 2_UCLA_Code/                       # Reference MATLAB implementation
│   └── ...                            #   Standard SSA tutorial (no gap support)
│
├── 3_Schoellhamer_Code/               # Adapted MATLAB implementation
│   ├── original/                      #   Schoellhamer's original code
│   └── modified/                      #   Our modifications for GNSS workflows
│
├── 4_Main_SSA_CERI/                   # Main Python pipeline
│   ├── SSA_Schoellhamer_UofM_V3.py    #   Core SSA implementation
│   ├── ToSSA_Input/                   #   Worked example with sample data
│   │   ├── input_TS/
│   │   └── NevGeoLab/
│   ├── SSA_OutPut/
│   │   ├── output_Figures/
│   │   └── output_TS/
│   └── README.txt                     #   Detailed usage guide
│
└── README.md                          # This file
```

---

## Installation

### Prerequisites
- Python **3.9+**
- Git

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/HadiHShali/SSA.git
cd SSA

# 2. (Recommended) Create and activate a virtual environment
python -m venv .venv

#    Windows (CMD):
.venv\Scripts\activate
#    macOS / Linux:
source .venv/bin/activate

# 3. Install dependencies
python -m pip install --upgrade pip
pip install numpy pandas requests scipy matplotlib hectorp
```

> The examples below assume the virtual environment is activated.

---

## Usage

### Option A — Run on your own (detrended) data

Use this path if you already have a GNSS time series with jumps and trend removed.

1. Place your input series in `4_Main_SSA_CERI/ToSSA_Input/input_TS/`.
2. Follow the step-by-step instructions in `4_Main_SSA_CERI/README.txt`.
3. Run the core pipeline:

```bash
cd 4_Main_SSA_CERI
python SSA_Schoellhamer_UofM_V3.py
```

Outputs (decomposed components, reconstructions, and figures) are written to `SSA_OutPut/output_TS/` and `SSA_OutPut/output_Figures/`.

### Option B — Start from Nevada Geodetic Laboratory data

This path pulls raw data straight from the [Nevada Geodetic Laboratory](http://geodesy.unr.edu/) and detrends it with Hectorp before running SSA.

**1. Download `.tenv3` files for a region of interest.**

The command below retrieves every station on the North American plate within `35° ≤ lat ≤ 36°`, `-90° ≤ lon ≤ -89°` that has at least 5 years of data, and saves the `.tenv3` files to a `tenv3/` folder:

```bash
python download_gnss_stations.py \
  --ref NA \
  --min-lat 35.0 --max-lat 36.0 \
  --min-lon -90 --max-lon -89 \
  --min-duration 5
```

**2. Parse the `.tenv3` files and build per-component directories.**

This reads the downloaded files and organizes each component (N/E/U) into its own directory under `Stns_Dir`:

```bash
python tenv3_reader_compDir.py
```

**3. Detrend, remove jumps, and run SSA.**

With the virtual environment activated, remove the trend and offsets with Hectorp, then run the SSA decomposition as described in `4_Main_SSA_CERI/README.txt`:

```bash
cd 4_Main_SSA_CERI
python SSA_Schoellhamer_UofM_V3.py
```

> *Note:* fill in the exact Hectorp detrending command you use here so others can reproduce the full pipeline end to end.

---

## References

- **Schoellhamer, D. H. (2001).** *Singular spectrum analysis for time series with missing data.* Geophysical Research Letters, 28(16), 3187–3190.
- **Hectorp** — [gitlab.com/machielsimonbos/hectorp](https://gitlab.com/machielsimonbos/hectorp)
- **Nevada Geodetic Laboratory** — [geodesy.unr.edu](http://geodesy.unr.edu/)

---

## Author

**Hadi H. Shali** — [@HadiHShali](https://github.com/HadiHShali)
<!-- Add your LinkedIn / email / website here to make it easy for collaborators and recruiters to reach you. -->

If you use this toolkit in your work, a citation or a link back is appreciated.

---

## License

<!-- Choose a license and add it here. MIT is a common, permissive default for research code. -->
*To be specified — add a `LICENSE` file and reference it here.*
