#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download GNSS station time series (.tenv3) data from UNR Geodesy Lab.

Workflow:
    1. If the reference frame file (midas.{Ref_Frame}.txt) is not in the current
       directory, download it automatically.
    2. Read the reference frame file and filter stations by lat/lon and duration.
    3. Download .tenv3 files for filtered stations in parallel.

Examples
--------
    # Default: North America plate, Memphis area, 3+ years of data
    python download_gnss_stations.py

    # Eurasia plate, Tehran area, 5+ years
    python download_gnss_stations.py --ref EU --min-lat 35.0 --max-lat 36.0 --min-lon 51.0 --max-lon 52.0 --min-duration 5

    # Pacific plate, custom output folder, more parallel threads
    python download_gnss_stations.py --ref PA --output ./data --workers 10

    # See all options
    python download_gnss_stations.py --help

Supported reference frames:
    IGS  - IGS20 (default global frame)
    AF   - Africa             NA - North America        SA - South America
    AN   - Antarctica         NB - North Bismark        SB - South Bismark
    AR   - Arabia             NZ - Nazca                SC - Scotia
    AU   - Australia          OK - Okhotsk              SL - Shetland
    BG   - Bering             ON - Okinawa              SO - Somalia
    BU   - Burma              PA - Pacific              SU - Sunda
    CA   - Caribbean          PM - Panama               WL - Woodlark
    CO   - Cocos              PS - Philippine Sea
    EU   - Eurasian           IN - Indian               MA - Mariana

Source: https://geodesy.unr.edu/gps_timeseries/IGS20/midas/

Created on Tue Jun 25 07:58:53 2024
Author: Hadi Heydarizadeh Shali
"""

import os
import sys
import time
import argparse
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed


# ─── VALID REFERENCE FRAME CODES ──────────────────────────────────────────────
VALID_REF_FRAMES = [
    "IGS", "AF", "AN", "AR", "AU", "BG", "BU", "CA", "CO", "EU",
    "IN", "MA", "NA", "NB", "NZ", "OK", "ON", "PA", "PM", "PS",
    "SA", "SB", "SC", "SL", "SO", "SU", "WL"
]


# ─── REFERENCE FRAME FILE DOWNLOADER ──────────────────────────────────────────
def download_reference_frame(ref_frame, output_dir):
    """
    Download the MIDAS reference frame file from UNR Geodesy Lab.

    If the file already exists locally, skip the download.

    Parameters
    ----------
    ref_frame : str
        Reference frame code (e.g., 'NA', 'EU', 'IGS').
    output_dir : str
        Directory where the file will be saved.

    Returns
    -------
    str
        Full path to the downloaded (or existing) reference frame file.
    """
    filename = f"midas.{ref_frame}.txt"
    output_path = os.path.join(output_dir, filename)
    url = f"https://geodesy.unr.edu/gps_timeseries/IGS20/midas/{filename}"

    if os.path.exists(output_path):
        print(f"[OK] Reference frame file already exists: {filename}")
        return output_path

    print(f"[..] Downloading reference frame file from:\n     {url}")
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"[OK] Saved to: {output_path}\n")
            return output_path
        else:
            raise RuntimeError(
                f"Failed to download reference frame file. "
                f"Status code: {response.status_code}"
            )
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Error downloading reference frame file: {e}")


# ─── STATION DATA DOWNLOADER ──────────────────────────────────────────────────
def download_station_data(station, ref_frame, output_dir):
    """
    Download a single GNSS station's .tenv3 time series file.

    If the file already exists locally, skip the download.

    Parameters
    ----------
    station : str
        4-character station ID (e.g., 'LCHS').
    ref_frame : str
        Reference frame code (e.g., 'NA').
    output_dir : str
        Directory where the file will be saved.
    """
    url = f"https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/{ref_frame}/{station}.{ref_frame}.tenv3"
    filename = f"{station}.{ref_frame}.tenv3.txt"
    output_path = os.path.join(output_dir, filename)

    if os.path.exists(output_path):
        print(f"[OK] {station}: already exists, skipping.")
        return

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"[..] {station}: downloaded -> {filename}")
        else:
            print(f"[XX] {station}: failed (status {response.status_code})")
    except requests.exceptions.RequestException as e:
        print(f"[XX] {station}: error - {e}")


# ─── PARALLEL DOWNLOAD ORCHESTRATOR ───────────────────────────────────────────
def parallel_download(stations, ref_frame, output_dir, max_workers):
    """
    Download multiple station files in parallel using a thread pool.
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_station = {
            executor.submit(download_station_data, station, ref_frame, output_dir): station
            for station in stations
        }
        for future in as_completed(future_to_station):
            station = future_to_station[future]
            try:
                future.result()
            except Exception as e:
                print(f"[XX] {station}: unexpected error - {e}")


# ─── COMMAND-LINE ARGUMENT PARSER ─────────────────────────────────────────────
def parse_arguments():
    """
    Parse command-line arguments. Returns an argparse.Namespace object.
    """
    parser = argparse.ArgumentParser(
        description="Download GNSS station time series (.tenv3) from UNR Geodesy Lab.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Default: North America plate, Memphis area
  python download_gnss_stations.py

  # Eurasia plate, Tehran area
  python download_gnss_stations.py --ref EU \\
      --min-lat 35.0 --max-lat 36.0 \\
      --min-lon 51.0 --max-lon 52.0

  # Pacific plate with custom output folder
  python download_gnss_stations.py --ref PA --output ./my_data
        """
    )

    parser.add_argument(
        "--ref",
        type=str,
        default="NA",
        choices=VALID_REF_FRAMES,
        metavar="CODE",
        help="Reference frame code (default: NA). Valid: " + ", ".join(VALID_REF_FRAMES)
    )

    parser.add_argument(
        "--min-lat", type=float, default=35.0,
        help="Minimum latitude in degrees (default: 35.0)"
    )
    parser.add_argument(
        "--max-lat", type=float, default=36.0,
        help="Maximum latitude in degrees (default: 36.0)"
    )
    parser.add_argument(
        "--min-lon", type=float, default=-90.0,
        help="Minimum longitude in degrees (default: -90.0)"
    )
    parser.add_argument(
        "--max-lon", type=float, default=-89.0,
        help="Maximum longitude in degrees (default: -89.0)"
    )

    parser.add_argument(
        "--min-duration", type=float, default=3.0,
        help="Minimum station observation duration in years (default: 3.0)"
    )

    parser.add_argument(
        "--output", type=str, default="tenv3",
        metavar="DIR",
        help="Output directory for downloaded .tenv3 files (default: ./tenv3)"
    )

    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of parallel download threads (default: 1)"
    )

    return parser.parse_args()


# ─── MAIN PIPELINE ────────────────────────────────────────────────────────────
def main():
    args = parse_arguments()
    t_start = time.time()

    # Validate inputs
    if args.min_lat >= args.max_lat:
        sys.exit("Error: --min-lat must be less than --max-lat.")
    if args.min_lon >= args.max_lon:
        sys.exit("Error: --min-lon must be less than --max-lon.")
    if args.min_duration < 0:
        sys.exit("Error: --min-duration must be non-negative.")
    if args.workers < 1:
        sys.exit("Error: --workers must be at least 1.")

    # 1) Set up directories
    base_dir = os.getcwd()
    output_directory = os.path.join(base_dir, args.output)
    os.makedirs(output_directory, exist_ok=True)

    # 2) Print summary of selected options
    print("=" * 70)
    print("  GNSS Station Downloader -- UNR Geodesy Lab")
    print("=" * 70)
    print(f"  Reference frame  : {args.ref}")
    print(f"  Lat range        : [{args.min_lat}, {args.max_lat}]")
    print(f"  Lon range        : [{args.min_lon}, {args.max_lon}]")
    print(f"  Min duration     : {args.min_duration} years")
    print(f"  Output directory : {output_directory}")
    print(f"  Parallel workers : {args.workers}")
    print("=" * 70 + "\n")

    # 3) Ensure the reference frame file is present (download if missing)
    ref_frame_path = download_reference_frame(args.ref, base_dir)

    # 4) Load and filter the reference frame file
    # The MIDAS file has NO header row — column names are defined manually
    # based on the official UNR Geodesy MIDAS README. The file has 27 columns:
    midas_columns = [
        'stnID',         # col 1  - 4-character station ID
        'version',       # col 2  - MIDAS version label
        'epoch_first',   # col 3  - time series first epoch (decimal year)
        'epoch_last',    # col 4  - time series last epoch (decimal year)
        'duration',      # col 5  - time series duration (years)
        'n_epochs_all',  # col 6  - number of epochs (used or not)
        'n_epochs_good', # col 7  - number of epochs of good data
        'n_pairs',       # col 8  - number of velocity sample pairs
        'Ve', 'Vn', 'Vu',          # cols 9-11  - east/north/up velocities (m/yr)
        'sVe', 'sVn', 'sVu',       # cols 12-14 - velocity uncertainties (m/yr)
        'Oe', 'On', 'Ou',          # cols 15-17 - offsets at first epoch (m)
        'fOe', 'fOn', 'fOu',       # cols 18-20 - fraction of outliers (e/n/u)
        'sdVe', 'sdVn', 'sdVu',    # cols 21-23 - std dev of velocity pairs
        'n_steps',       # col 24 - assumed number of steps
        'Lat',           # col 25 - station latitude (degrees)
        'Lon',           # col 26 - station longitude (degrees)
        'Hgt',           # col 27 - station height (meters)
    ]

    df = pd.read_csv(
        ref_frame_path,
        sep=r'\s+',
        header=None,
        names=midas_columns,
        comment='#',
    )

    filtered_df = df[
        (df['Lat'] >= args.min_lat) & (df['Lat'] <= args.max_lat) &
        (df['Lon'] >= args.min_lon) & (df['Lon'] <= args.max_lon) &
        (df['duration'] >= args.min_duration)
    ]
    resultant_stations = filtered_df['stnID'].unique().tolist()

    if not resultant_stations:
        print("[!] No stations match the filter criteria. Exiting.")
        return

    print(f"[i] Found {len(resultant_stations)} stations matching criteria:")
    print(f"    {', '.join(resultant_stations)}\n")

    # 5) Download station data in parallel
    parallel_download(
        resultant_stations,
        args.ref,
        output_directory,
        max_workers=args.workers
    )

    # 6) Summary
    t_end = time.time()
    print(f"\n[OK] Done. Downloaded data for {len(resultant_stations)} stations "
          f"in {t_end - t_start:.1f} seconds.")


if __name__ == "__main__":
    main()