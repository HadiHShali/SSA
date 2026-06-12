# -*- coding: utf-8 -*-
"""
Singular Spectrum Analysis (SSA) for GNSS time series.

Example
-------
    python SSA_Schoellhamer_UofM_V3.py --MaxEig 16 --plot-reconstruct --plot-EigValue --plot-lomb

    python SSA_Schoellhamer_UofM_V3.py \\
        --input  ToSSA_Input/input_TS/ \\
        --output-ts  SSA_OutPut/output_TS/ \\
        --output-fig SSA_OutPut/output_Figures/ \\
        --embedding 910 --MaxEig 16 \\
        --components "1-2,3-4,5-6" \\
        --plot-all
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz
from astropy.timeseries import LombScargle
import os


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="SSA processing for GNSS time series",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # I/O
    p.add_argument("--input", default="ToSSA_Input/input_TS/",
                   help="Directory containing .dat input files")
    p.add_argument("--output-ts", dest="output_ts",
                   default="SSA_OutPut/output_TS/",
                   help="Directory for SSA cyclic time-series outputs (.mom)")
    p.add_argument("--output-fig", dest="output_fig",
                   default="SSA_OutPut/output_Figures/",
                   help="Root figure directory (one sub-folder per station)")

    # Processing
    p.add_argument("--embedding", "-M", type=int, default=910, metavar="M",
                   help="Embedding dimension (days); ~2.5 yr = 910")
    p.add_argument("--MaxEig", type=int, default=16, metavar="N",
                   help="Number of eigenvalues used in reconstruction")
    p.add_argument("--components", default="1-2,3-4,5-6",
                   help="RC groupings for multi-panel plot, e.g. '1-2,3-4,5-6'")
    p.add_argument("--toeplitz", action="store_true",
                   help="Use Toeplitz covariance instead of trajectory-matrix covariance")
    p.add_argument("--stations", nargs="+", default=None, metavar="STA",
                   help="Process only these station names (default: all in --input)")
    p.add_argument("--dry-run", action="store_true",
                   help="List stations that would be processed, then exit")

    # Lomb-Scargle
    p.add_argument("--min-freq", dest="min_freq", type=float, default=0.3,
                   help="Lomb-Scargle lower bound (cycles / year)")
    p.add_argument("--max-freq", dest="max_freq", type=float, default=100.0,
                   help="Lomb-Scargle upper bound (cycles / year)")

    # Plot flags
    p.add_argument("--plot-all", dest="plot_all", action="store_true",
                   help="Enable every available plot")
    p.add_argument("--plot-gappy", dest="plot_gappy", action="store_true",
                   help="Gappy input time series (mean-removed, normalised)")
    p.add_argument("--plot-toeplitz", dest="plot_toeplitz", action="store_true",
                   help="Toeplitz diagonal values and lagged-correlation matrix")
    p.add_argument("--plot-EigValue", dest="plot_eigvalue", action="store_true",
                   help="Eigenvalue spectrum — all in blue, selected in red")
    p.add_argument("--plot-reconstruct", dest="plot_reconstruct", action="store_true",
                   help="Reconstructed signal (sum of RC 1–MaxEig) overlaid on data")
    p.add_argument("--plot-individuals", dest="plot_individuals", action="store_true",
                   help="One figure per RC component from 1 to MaxEig")
    p.add_argument("--plot-residuals", dest="plot_residuals", action="store_true",
                   help="Residuals (data − reconstruction)")
    p.add_argument("--plot-lomb", dest="plot_lomb", action="store_true",
                   help="Lomb-Scargle PSD of reconstruction and residuals")

    return p.parse_args()


# =============================================================================
# Helpers
# =============================================================================
def resolve_plots(args):
    """Propagate --plot-all to every individual flag."""
    if args.plot_all:
        for attr in ("plot_gappy", "plot_toeplitz", "plot_eigvalue",
                     "plot_reconstruct", "plot_individuals",
                     "plot_residuals", "plot_lomb"):
            setattr(args, attr, True)
    return args


def parse_components(comp_str, maxeig):
    """
    '1-2,3-6,7-15' → [(0,2),(2,6),(6,15)]  (0-based start, exclusive end).
    Each group's end is clipped to maxeig so you never index beyond what
    was reconstructed, but the label shown in the plot reflects the numbers
    the user typed, not the clipped values.

    Returns list of (start, end, label) tuples.
    """
    groups = []
    for part in comp_str.split(","):
        a, b = part.strip().split("-")
        a, b = int(a), int(b)
        start = a - 1                   # 0-based
        end   = min(b, maxeig)          # clip to available RCs
        label = f"RC {a}–{end}"        # label reflects actual range used
        if start >= maxeig:
            continue                    # entirely out of range — skip silently
        groups.append((start, end, label))
    return groups


def savefig(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_psd_LS(raw_values, fmin_cpy=0.3, fmax_cpy=100.0):
    """
    Lomb-Scargle PSD for a daily time series with gaps (NaNs).
    Unlike Welch, Lomb-Scargle operates on the valid (time, value) pairs
    directly and handles non-uniform sampling / missing data natively.

    Parameters
    ----------
    raw_values : 1-D array  (daily samples; NaNs for gaps)
    fmin_cpy, fmax_cpy : frequency window in cycles per year

    Returns
    -------
    freq_cpy : ndarray, frequencies in cycles per year
    psd      : ndarray, Lomb-Scargle PSD (normalisation="psd")
    """
    x = np.asarray(raw_values, dtype=float)
    mask = ~np.isnan(x)
    t_days = np.arange(len(x), dtype=float)[mask]
    vals = x[mask] - np.nanmean(x[mask])

    fmin_cpd = fmin_cpy / 365.25
    fmax_cpd = fmax_cpy / 365.25

    ls = LombScargle(t_days, vals, normalization="psd")
    freq_cpd, psd = ls.autopower(
        minimum_frequency=fmin_cpd,
        maximum_frequency=fmax_cpd,
        samples_per_peak=5,
    )
    freq_cpy = freq_cpd * 365.25
    return freq_cpy, psd


# =============================================================================
# SSA core
# =============================================================================
def SSABasicCleanUCLAtoeplitz(trajmat, toep=None, maxieg=None):
    """
    Perform Singular Spectrum Analysis (SSA) to clean the trajectory matrix
    with or without Toeplitz covariance.

    Parameters
    ----------
    trajmat : ndarray  (L x M)  trajectory matrix
    toep    : optional Toeplitz covariance matrix
    maxieg  : optional cap on number of eigenvectors used

    Returns
    -------
    dict with keys RC, LAMBDA, RHO, num_nans
    """
    print(f"Entering SSABasicCleanUCLAtoeplitz, inputs: {len(locals())}")

    L, M = trajmat.shape
    originalvec = np.concatenate((trajmat[:, 0], trajmat[-1, 1:]))
    num_nans = np.sum(np.isnan(originalvec))
    N = L + M - 1

    if num_nans != 0:
        print("NaNs detected in the input matrix.")
        trajmat_ele_good = np.copy(trajmat)
        trajmat_ele_good[np.isnan(trajmat_ele_good)] = 0
        trajmat_ele_TF = np.ones_like(trajmat)
        trajmat_ele_TF[np.isnan(trajmat)] = 0

        pct_good = np.sum(trajmat_ele_TF) / trajmat.size
        print(f"Percentage of good data: {pct_good:.4f}")

        C = trajmat_ele_good.T @ trajmat_ele_good
        NaNadj = trajmat_ele_TF.T @ trajmat_ele_TF
        C = C / NaNadj

        if np.any(np.isinf(C)):
            print("Bad elements in covariance matrix (inf / divide-by-zero)")
            return
    else:
        print("No NaNs detected, using normal covariance")
        trajmat_ele_good = trajmat
        C = trajmat.T @ trajmat

    if toep is not None:
        C = toep
        print("Using Toeplitz matrix for covariance")

    try:
        LAMBDAM, RHO = np.linalg.eig(C)
    except np.linalg.LinAlgError:
        print("Eigenvalue/eigenvector calculation failed.")
        return

    LAMBDA = LAMBDAM.copy()
    LAMBDA[LAMBDA < 0] = 1e-16

    ind = np.argsort(LAMBDA)[::-1]
    LAMBDA = LAMBDA[ind]
    RHO = RHO[:, ind]

    if maxieg is not None:
        maxieg = min(maxieg, len(LAMBDA))
        RHO_use = RHO[:, :maxieg]
        PC = trajmat_ele_good @ RHO_use
        if np.any(np.isnan(PC)):
            print("NaNs found in PC — fix required")
            return
        print("No NaNs in PCs, continuing with reconstruction.")
        RC = np.zeros((N, maxieg))
        for m in range(maxieg):
            buf = np.outer(PC[:, m], RHO_use[:, m])
            buf = np.flipud(buf)
            for n in range(N):
                RC[n, m] = np.mean(np.diagonal(buf, offset=-(N - M) + n))
    else:
        PC = trajmat_ele_good @ RHO
        if np.any(np.isnan(PC)):
            print("NaNs found in PC — fix required")
            return
        print("No NaNs in PCs, continuing with reconstruction.")
        RC = np.zeros((N, M))
        for m in range(M):
            buf = np.outer(PC[:, m], RHO[:, m])
            buf = np.flipud(buf)
            for n in range(N):
                RC[n, m] = np.mean(np.diagonal(buf, offset=-(N - M) + n))

    print("Leaving SSABasicCleanUCLAtoeplitz")
    return {"RC": RC, "LAMBDA": LAMBDA, "RHO": RHO, "num_nans": num_nans}


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    args = resolve_plots(args)

    input_dir  = args.input
    output_ts  = args.output_ts
    output_fig = args.output_fig
    M          = args.embedding
    maxeig     = args.MaxEig

    os.makedirs(output_ts,  exist_ok=True)
    os.makedirs(output_fig, exist_ok=True)

    # Collect files to process
    file_names = [f for f in os.listdir(input_dir)
                  if os.path.isfile(os.path.join(input_dir, f))]

    if args.stations:
        file_names = [f for f in file_names
                      if any(s in f for s in args.stations)]

    if args.dry_run:
        print("Stations that would be processed:")
        for f in file_names:
            print(f"  {f}")
        return

    comp_groups = parse_components(args.components, maxeig)

    # ------------------------------------------------------------------
    for fldrs in file_names:
        St_Names = fldrs[:10]
        Compnnt  = fldrs[5]
        Comp = {"0": "E", "1": "N", "2": "U"}.get(Compnnt, "?")

        print(f"\n{'='*60}")
        print(f"Processing: {St_Names}  component: {Comp}")

        # Per-station figure folder
        fig_dir = os.path.join(output_fig, St_Names)
        os.makedirs(fig_dir, exist_ok=True)

        # ---- Load data ---------------------------------------------------
        filename_dir = os.path.join(input_dir, f"{St_Names}.dat")
        comments, data_lines = [], []
        with open(filename_dir) as fh:
            for line in fh:
                if line.startswith("#"):
                    comments.append(line.strip())
                else:
                    data_lines.append(line)

        full_data = np.loadtxt(data_lines)
        n_data    = full_data.shape[0]

        tmjd   = full_data[:, 0]   # Modified Julian Date
        tyfrac = full_data[:, 1]   # Fractional year
        obser  = full_data[:, 2]   # Raw observation (mm)
        x_in   = full_data[:, 3]   # Detrended, jumps removed (mm)

        dataindx = tmjd - tmjd[0] + 1

        x_std  = np.std(x_in, ddof=1)
        x_mean = np.mean(x_in)
        x_norm = (x_in - x_mean) / x_std

        # Continuous time vector
        t = np.arange(tmjd[0], tmjd[-1] + 1) - tmjd[0] + 1
        N = len(t)

        if N == n_data:
            print("No gaps in input data")
        else:
            print(f"Gaps: {N - n_data} missing days")

        if n_data < 2 * M:
            print("Time series shorter than 2*M — skipping.")
            continue

        # Build gappy arrays
        idx = dataindx.astype(int) - 1

        x_in_NaN   = np.full(N, np.nan)
        x_in_NaN[idx] = x_in

        x_norm_NaN = np.full(N, np.nan)
        x_norm_NaN[idx] = x_norm

        tFrac_NaN  = np.full(N, np.nan)
        tFrac_NaN[idx] = tyfrac

        X = np.where(np.isnan(x_norm_NaN), 0.0, x_norm_NaN)  # NaNs → 0

        NaNindx  = np.where(np.isnan(x_norm_NaN))[0]
        num_nans = len(NaNindx)
        print(f"Data pts: {n_data}, continuous pts: {N}, "
              f"gaps: {num_nans} ({100*num_nans/N:.1f} %)")

        # ---- Trajectory matrix -------------------------------------------
        l        = N - M + 1
        rowind   = np.arange(l)[:, None]
        colind   = np.arange(1, M + 1)
        trajmat  = x_norm_NaN[colind + rowind - 1]

        # ---- Toeplitz covariance ------------------------------------------
        Cloop = np.zeros(M)
        for J in range(M):
            nterms = N - J
            products = X[:nterms] * X[J:J + nterms]
            nonzero  = products[products != 0]
            Cloop[J] = nonzero.mean() if len(nonzero) else 0.0

        Ctoep = toeplitz(Cloop)

        # ---- Plot: gappy input -------------------------------------------
        if args.plot_gappy:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.plot(tFrac_NaN, x_in_NaN, "bo-", ms=2, lw=0.3)
            ax.set_xlabel("Time (fractional year)")
            ax.set_ylabel("Position (mm)")
            ax.set_title(f"{St_Names} — gappy input (detrended, jumps removed)")
            savefig(fig, os.path.join(fig_dir, "01_gappy_input.png"))

        # ---- Plot: Toeplitz ----------------------------------------------
        if args.plot_toeplitz:
            fig, axes = plt.subplots(1, 2, figsize=(14, 5))
            axes[0].plot(Cloop, "b-")
            axes[0].set_xlabel("Lag (days)")
            axes[0].set_ylabel("Correlation")
            axes[0].set_title(f"{St_Names} — Toeplitz diagonal")
            im = axes[1].imshow(Ctoep, aspect="auto", cmap="viridis")
            fig.colorbar(im, ax=axes[1], label="Correlation")
            axes[1].set_title(f"{St_Names} — Toeplitz matrix  (M={M})")
            savefig(fig, os.path.join(fig_dir, "02_toeplitz.png"))

        # ---- Run SSA -----------------------------------------------------
        toep_arg = Ctoep if args.toeplitz else None
        out = SSABasicCleanUCLAtoeplitz(trajmat, toep=toep_arg, maxieg=maxeig)
        if out is None:
            print("SSA failed — skipping station.")
            continue

        RC = out["RC"]
        LAMBDA = out["LAMBDA"]

        # Re-insert NaNs at gap positions
        RCNaN = RC.copy()
        RCNaN[NaNindx, :] = np.nan

        # Reconstructed signal (sum of first maxeig RCs), back in mm
        RC_sum     = np.sum(RCNaN[:, :maxeig], axis=1)
        RC_sum_mm  = RC_sum * x_std + x_mean   # unnormalised

        # Residual (mm)
        residual_mm = x_in_NaN - RC_sum_mm

        # ---- Plot: eigenvalue spectrum -----------------------------------
        if args.plot_eigvalue:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.semilogy(LAMBDA, "b+-", label="all eigenvalues")
            ax.semilogy(np.arange(maxeig), LAMBDA[:maxeig], "ro-",
                        label=f"selected (1–{maxeig})")
            ax.set_xlabel("Index")
            ax.set_ylabel("Eigenvalue")
            ax.set_title(f"{St_Names} — eigenvalue spectrum")
            ax.legend()
            savefig(fig, os.path.join(fig_dir, "03_eigenvalues.png"))

        # ---- Plot: reconstruction vs data --------------------------------
        if args.plot_reconstruct:
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.plot(tFrac_NaN, x_in_NaN,  lw=0.8, label="data")
            ax.plot(tFrac_NaN, RC_sum_mm, "--", lw=1.2,
                    label=f"reconstruction (RC 1–{maxeig})")
            ax.set_xlabel("Time (fractional year)")
            ax.set_ylabel("Position (mm)")
            ax.set_title(f"{St_Names} — reconstruction vs data")
            ax.legend()
            savefig(fig, os.path.join(fig_dir, "04_reconstruction.png"))

        # ---- Plot: multi-group panel -------------------------------------
        if args.plot_reconstruct and comp_groups:
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.plot(tFrac_NaN, x_norm_NaN, lw=0.6, label="data (norm.)")
            for (a, b, label) in comp_groups:
                ax.plot(tFrac_NaN,
                        np.sum(RCNaN[:, a:b], axis=1),
                        lw=1.0, label=label)
            ax.set_xlabel("Time (fractional year)")
            ax.set_ylabel("Normalised amplitude")
            ax.set_title(f"{St_Names} — RC groups  (M={M}, MaxEig={maxeig})")
            ax.legend()
            savefig(fig, os.path.join(fig_dir, "05_rc_groups.png"))

        # ---- Plot: individual RC components ------------------------------
        if args.plot_individuals:
            for k in range(min(maxeig, RC.shape[1])):
                fig, ax = plt.subplots(figsize=(12, 3))
                ax.plot(tFrac_NaN, RCNaN[:, k] * x_std + x_mean,
                        lw=0.8, color="steelblue")
                ax.set_xlabel("Time (fractional year)")
                ax.set_ylabel("Position (mm)")
                ax.set_title(f"{St_Names} — RC {k+1}")
                savefig(fig, os.path.join(fig_dir,
                        f"06_rc_{k+1:02d}.png"))

        # ---- Plot: residuals ---------------------------------------------
        if args.plot_residuals:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.plot(tFrac_NaN, residual_mm, lw=0.6, color="gray")
            ax.axhline(0, color="k", lw=0.5, ls="--")
            ax.set_xlabel("Time (fractional year)")
            ax.set_ylabel("Residual (mm)")
            ax.set_title(f"{St_Names} — residuals (data − reconstruction)")
            savefig(fig, os.path.join(fig_dir, "07_residuals.png"))

        # ---- Plot: Lomb-Scargle PSD --------------------------------------
        if args.plot_lomb:
            freq_rc,  psd_rc  = compute_psd_LS(RC_sum_mm,
                                                fmin_cpy=args.min_freq,
                                                fmax_cpy=args.max_freq)
            freq_res, psd_res = compute_psd_LS(residual_mm,
                                                fmin_cpy=args.min_freq,
                                                fmax_cpy=args.max_freq)
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.loglog(freq_rc,  psd_rc,  "r",   lw=0.8, label="cyclics (RC sum)")
            ax.loglog(freq_res, psd_res, "b.-", ms=1,   label="residuals / noise")
            ax.set_xlabel("Frequency (cycles / year)")
            ax.set_ylabel("PSD")
            ax.set_title(f"{St_Names} — Lomb-Scargle PSD  "
                         f"[{args.min_freq}–{args.max_freq} cpy]")
            ax.legend()
            ax.grid(True, which="both", ls=":")
            savefig(fig, os.path.join(fig_dir, "08_lomb_scargle.png"))

        # ---- Save SSA cyclic time series ---------------------------------
        os.makedirs(output_ts, exist_ok=True)
        valid = ~np.isnan(RC_sum_mm)
        out_ts_path = os.path.join(output_ts, f"{St_Names}_Cyclics.mom")
        with open(out_ts_path, "w", newline="") as fh:
            for ti, ri in zip(tmjd, RC_sum_mm[valid]):
                fh.write(f"{ti} {ri:.6f}\n")

        # ---- Save residual .mom for Hector -------------------------------
        # directory = os.getcwd()
        # output_dir_mom = os.path.join(directory, "..",
                                      # "4th_Velocity_HectorP_SSACycleRmved")
        # subfolder_mom = os.path.join(output_dir_mom, "obs_files")
        # os.makedirs(subfolder_mom, exist_ok=True)

        # sampling_period = "# sampling period 1.0"
        # out_mom_path = os.path.join(subfolder_mom, f"{St_Names}_{Comp}.mom")
        # RC_filt = RC_sum_mm[valid]
        # with open(out_mom_path, "w", newline="") as fh:
            # fh.write(sampling_period + "\n")
            # if comments and comments[0].startswith("# exp"):
                # comments[0] = comments[0].replace("# exp", "# log", 1)
            # for comment in comments:
                # fh.write(comment + "\n")
            # for ti, oi, ri in zip(tmjd, obser, RC_filt):
                # fh.write(f"{ti} {oi - ri:.6f}\n")

        print(f"  -> TS saved : {out_ts_path}")
        #print(f"  -> MOM saved: {out_mom_path}")
        print(f"  -> Figures  : {fig_dir}/")

    print("\nDone.")


if __name__ == "__main__":
    main()
