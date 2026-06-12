#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tenv3_to_mom.py
---------------
Converts Nevada Geodetic Lab .tenv3 GNSS time series files into
Hector-compatible .mom observation files (East, North, Up components).

Reads jump/offset information from steps.txt and writes per-component
.mom files with proper offset and log headers.

Author: Hadi Heydarizadeh Shali
"""

import os
import datetime
import urllib.request
import pandas as pd


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TENV3_COLUMNS = [
    'site', 'YYMMMDD', 'yyyy.yyyy', '__MJD', 'week', 'd',
    'reflon', '_e0(m)', '__east(m)', '____n0(m)', '_north(m)',
    'u0(m)', '____up(m)', '_ant(m)', 'sig_e(m)', 'sig_n(m)',
    'sig_u(m)', '__corr_en', '__corr_eu', '__corr_nu',
    '_latitude(deg)', '_longitude(deg)', '__height(m)'
]

MJD_EPOCH_OFFSET = 678576  # Ordinal offset for Modified Julian Date


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def ensure_dir(path: str) -> None:
    """Create directory if it does not exist."""
    os.makedirs(path, exist_ok=True)


STEPS_URL = 'https://geodesy.unr.edu/NGLStationPages/steps.txt'


def ensure_steps_file(steps_file: str) -> None:
    """
    Check if steps.txt exists at the given path.
    If not, download it from the NGL server and save it locally.
    """
    if os.path.isfile(steps_file):
        print(f"[steps] Found '{steps_file}'.")
        return

    print(f"[steps] '{steps_file}' not found. Downloading from:\n  {STEPS_URL}")
    try:
        urllib.request.urlretrieve(STEPS_URL, steps_file)
        size_kb = os.path.getsize(steps_file) / 1024
        print(f"[steps] Downloaded successfully ({size_kb:.0f} KB).")
    except Exception as e:
        raise RuntimeError(
            f"Could not download steps.txt from {STEPS_URL}.\n"
            f"Please download it manually and place it at '{steps_file}'.\n"
            f"Error: {e}"
        )


def date_str_to_mjd(date_str: str) -> int:
    """Convert a date string (format: YYMONDD, e.g. '10JAN15') to MJD."""
    date_obj = datetime.datetime.strptime(date_str, '%y%b%d')
    return date_obj.toordinal() - MJD_EPOCH_OFFSET


def load_steps(steps_file: str) -> list[list[str]]:
    """Parse the steps.txt file into a list of tokenized rows."""
    with open(steps_file, 'r') as f:
        return [line.split() for line in f.readlines()]


def get_station_steps(steps: list, station: str) -> list[list]:
    """
    Filter steps for a given station and attach the MJD column.
    Returns a list of rows with MJD inserted at index 2.
    """
    station_steps = [row[:] for row in steps if row[0] == station]
    for row in station_steps:
        row.insert(2, date_str_to_mjd(row[1]))
    return station_steps


def read_tenv3(filepath: str) -> pd.DataFrame | None:
    """
    Read a .tenv3 file into a DataFrame.
    Returns None if the file cannot be decoded.
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()[1:]  # skip header
        data = [line.split() for line in lines]
        df = pd.DataFrame(data, columns=TENV3_COLUMNS)
        numeric_cols = df.columns[df.columns.get_loc('YYMMMDD') + 1:]
        df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')
        return df
    except UnicodeDecodeError as e:
        print(f"[WARNING] Cannot decode {filepath}: {e}")
        return None


def compute_displacement(df: pd.DataFrame) -> tuple:
    """
    Compute mean-centred East, North, Up displacements in mm
    from tenv3 reference + residual columns.
    Returns (MJD series, East mm, North mm, Up mm).
    """
    mjd   = df['__MJD']
    east  = (df['_e0(m)']  + df['__east(m)']  - (df['_e0(m)']  + df['__east(m)'] ).mean()) * 1000
    north = (df['____n0(m)'] + df['_north(m)'] - (df['____n0(m)'] + df['_north(m)']).mean()) * 1000
    up    = (df['u0(m)']   + df['____up(m)']  - (df['u0(m)']   + df['____up(m)'] ).mean()) * 1000
    return mjd, east, north, up


def classify_offsets(station_steps: list, mjd_values) -> tuple[list, list]:
    """
    Separate steps into equipment-change offsets (flag=1) and
    seismic offsets with log decay (flag=2).

    If a step MJD is not in the data, it is advanced day-by-day
    until a matching epoch is found.

    Returns (eqc_offsets, eq_offsets).
    """
    mjd_set = set(mjd_values)
    eqc_offsets, eq_offsets = [], []

    for row in station_steps:
        mjd  = row[2]
        flag = int(row[3])

        # Advance MJD until it lands on an observed epoch
        while mjd not in mjd_set:
            mjd += 1

        target = eqc_offsets if flag == 1 else eq_offsets
        if mjd not in target:
            target.append(mjd)

    return eqc_offsets, eq_offsets


# ---------------------------------------------------------------------------
# File writing
# ---------------------------------------------------------------------------

def write_mom_header(f, eqc_offsets: list, eq_offsets: list) -> None:
    """Write the Hector .mom header block to an open file handle."""
    f.write('# sampling period 1.0\n')
    for mjd in eqc_offsets:
        f.write(f'# offset {mjd}\n')
    for mjd in eq_offsets:
        f.write(f'# offset {mjd}\n')
    for mjd in eq_offsets:
        f.write(f'# log  {mjd}  50\n')


def write_mom_files(output_folder: str, station: str,
                    mjd, east, north, up,
                    eqc_offsets: list, eq_offsets: list) -> None:
    """
    Write East (_0.mom), North (_1.mom), and Up (_2.mom) observation files
    for the given station under output_folder/<STATION>_<C>/obs_files/.
    """
    components = {
        'E': (0, east),
        'N': (1, north),
        'U': (2, up),
    }

    for comp, (idx, data) in components.items():
        out_dir = os.path.join(output_folder, f'{station}_{comp}', 'obs_files')
        ensure_dir(out_dir)
        filepath = os.path.join(out_dir, f'{station}_{idx}.mom')

        with open(filepath, 'w') as f:
            write_mom_header(f, eqc_offsets, eq_offsets)
            for t, val in zip(mjd, data):
                f.write(f'{t} {val:.6f}\n')


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def process_all_stations(tenv3_dir: str, output_dir: str,
                         steps_file: str = 'steps.txt') -> None:
    """
    Process every .tenv3 file in tenv3_dir and write .mom files
    into the output_dir hierarchy.
    """
    ensure_steps_file(steps_file)
    steps = load_steps(steps_file)
    ensure_dir(output_dir)

    station_files = os.listdir(tenv3_dir)
    print(f"Found {len(station_files)} station file(s) in '{tenv3_dir}'.\n")

    for filename in station_files:
        station = filename[:4]
        filepath = os.path.join(tenv3_dir, filename)

        df = read_tenv3(filepath)
        if df is None:
            continue

        station_steps          = get_station_steps(steps, station)
        mjd, east, north, up   = compute_displacement(df)
        eqc_offsets, eq_offsets = classify_offsets(station_steps, mjd.values)

        write_mom_files(output_dir, station, mjd, east, north, up,
                        eqc_offsets, eq_offsets)

        print(f"  [{station}] written  "
              f"(eqc (equipment changes)={len(eqc_offsets)}, eq (Earthquakes)={len(eq_offsets)})")

    print("\nDone.")


if __name__ == '__main__':
    BASE_DIR      = os.getcwd()
    TENV3_DIR     = os.path.join(BASE_DIR, 'tenv3')
    OUTPUT_DIR    = os.path.join(BASE_DIR, 'Stns_Dir')
    STEPS_FILE    = os.path.join(BASE_DIR, 'steps.txt')

    process_all_stations(TENV3_DIR, OUTPUT_DIR, STEPS_FILE)