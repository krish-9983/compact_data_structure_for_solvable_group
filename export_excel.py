#!/usr/bin/env python3
"""
export_excel.py — Build the three Excel files needed by plot_6.py.

Reads the raw benchmark CSVs and produces:
  <out>_plot1_FixedL.xlsx        — one row per group, fixed L (smallest L)
  <out>_plot2_PerLSize.xlsx      — one row per (group, L) combination
  <out>_plot3_GrandAverage.xlsx  — one row per group, averaged over all L values

Usage:
  python export_excel.py
  python export_excel.py --gap results/benchmark_results.csv \
                         --cpp results/benchmark_results_c++.csv \
                         --out results/benchmark_results
"""

import argparse
import re
import warnings
from pathlib import Path

import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")


GAP_COLS = [
    'Engine', 'Group', 'Order', 'Instructions',
    'Seed', 'Repetitions', 'Result', 'TotalTime_ms', 'TimePerRep_ms'
]


def fix_gap_name(name):
    """Convert GAP group name to standard format."""
    m = re.match(r'SmallGroup(\d+)_(\d+)$', name)
    if m:
        return f"SmallGroup({m.group(1)},{m.group(2)})"
    m = re.match(r'DirectProduct_?SmallGroup(\d+)_(\d+)_?CyclicGroup(\d+)', name)
    if m:
        return (f"DirectProduct(SmallGroup({m.group(1)},{m.group(2)}),"
                f"CyclicGroup({m.group(3)}))")
    return name


def short_name(name):
    """SmallGroup(64,267) -> 64x267"""
    m = re.search(r'\((\d+),(\d+)\)', name)
    return f"{m.group(1)}x{m.group(2)}" if m else name


def family(name):
    n = name.lower()
    if "directproduct" in n: return "DirectProduct"
    if "small"         in n: return "SmallGroup"
    return "Other"


def load(gap_csv, cpp_csv):
    """Load and merge GAP and CDS CSVs. Returns merged DataFrame."""
    cpp = pd.read_csv(cpp_csv)
    cpp.columns = cpp.columns.str.strip()

    # GAP CSV may or may not have a header
    with open(gap_csv) as f:
        first = f.readline()
    if first.lower().startswith("engine"):
        gap = pd.read_csv(gap_csv)
    else:
        gap = pd.read_csv(gap_csv, header=None, names=GAP_COLS)

    gap['Group'] = gap['Group'].apply(fix_gap_name)

    # Per-operation time in ms
    cpp['tpi'] = cpp['TimePerRep_ms'] / cpp['Instructions']
    gap['tpi'] = gap['TimePerRep_ms'] / gap['Instructions']

    # Rename L column for consistency
    cpp = cpp.rename(columns={'Instructions': 'L'})
    gap = gap.rename(columns={'Instructions': 'L'})

    return cpp, gap


def write_header(ws, title):
    """Write a bold title in row 1; data starts at row 3 (skiprows=2)."""
    ws.write(0, 0, title)


# ── Plot 1: fixed L (smallest L), one row per group ──────────────────────────

def make_plot1(cpp, gap, out_path):
    L_fixed = int(cpp['L'].min())
    n_seeds = cpp.groupby('Group')['Seed'].nunique().max()

    cpp1 = cpp[cpp['L'] == L_fixed]
    gap1 = gap[gap['L'] == L_fixed]

    cpp_avg = (cpp1.groupby('Group')['tpi']
                   .agg(['mean', 'std', lambda x: x.std()/x.mean()*100
                         if x.mean() > 0 else 0, 'min', 'max'])
                   .reset_index())
    cpp_avg.columns = ['Group', 'CDS avg time/instr (ms)', 'CDS std (ms)',
                       'CDS CV (%)', 'CDS min (ms)', 'CDS max (ms)']

    gap_avg = (gap1.groupby('Group')['tpi']
                   .agg(['mean', 'std', lambda x: x.std()/x.mean()*100
                         if x.mean() > 0 else 0, 'min', 'max'])
                   .reset_index())
    gap_avg.columns = ['Group', 'GAP avg time/instr (ms)', 'GAP std (ms)',
                       'GAP CV (%)', 'GAP min (ms)', 'GAP max (ms)']

    df = pd.merge(cpp_avg, gap_avg, on='Group')

    # Add order and family from CDS data
    meta = cpp1.groupby('Group').agg(Order=('Order', 'first')).reset_index()
    df   = pd.merge(df, meta, on='Group')
    df['Group (short)'] = df['Group'].apply(short_name)
    df['Group (full)']  = df['Group']
    df['Family']        = df['Group'].apply(family)
    df['L (fixed)']     = L_fixed
    df['Seeds used']    = n_seeds
    df['Speedup (GAP/CDS)'] = (
        df['GAP avg time/instr (ms)'] / df['CDS avg time/instr (ms)']
    )
    df = df.sort_values('Order')

    cols = [
        'Group (short)', 'Group (full)', 'Order', 'Family',
        'L (fixed)', 'Seeds used',
        'CDS avg time/instr (ms)', 'CDS std (ms)', 'CDS CV (%)',
        'CDS min (ms)', 'CDS max (ms)',
        'GAP avg time/instr (ms)', 'GAP std (ms)', 'GAP CV (%)',
        'GAP min (ms)', 'GAP max (ms)',
        'Speedup (GAP/CDS)'
    ]

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        pd.DataFrame([['Fixed-L benchmark results']]).to_excel(
            writer, sheet_name='Plot1_FixedL', index=False, header=False)
        df[cols].to_excel(
            writer, sheet_name='Plot1_FixedL',
            startrow=2, index=False)

    print(f"  saved  {Path(out_path).name}  ({len(df)} groups, L={L_fixed})")


# ── Plot 2: per (group, L), one row per combination ──────────────────────────

def make_plot2(cpp, gap, out_path):
    cpp_avg = (cpp.groupby(['Group', 'Order', 'L'])['tpi']
                  .agg(['mean', 'std', lambda x: x.std()/x.mean()*100
                        if x.mean() > 0 else 0])
                  .reset_index())
    cpp_avg.columns = ['Group', 'Order', 'L',
                       'CDS avg time/instr (ms)', 'CDS std (ms)', 'CDS CV (%)']

    gap_avg = (gap.groupby(['Group', 'L'])['tpi']
                  .agg(['mean', 'std', lambda x: x.std()/x.mean()*100
                        if x.mean() > 0 else 0])
                  .reset_index())
    gap_avg.columns = ['Group', 'L',
                       'GAP avg time/instr (ms)', 'GAP std (ms)', 'GAP CV (%)']

    df = pd.merge(cpp_avg, gap_avg, on=['Group', 'L'])
    df['Group (short)'] = df['Group'].apply(short_name)
    df['Group (full)']  = df['Group']
    df['Family']        = df['Group'].apply(family)
    df['Seeds used']    = cpp.groupby(['Group', 'L'])['Seed'].nunique().max()
    df['Speedup (GAP/CDS)'] = (
        df['GAP avg time/instr (ms)'] / df['CDS avg time/instr (ms)']
    )
    df = df.rename(columns={'L': 'L (SLP size)'})
    df = df.sort_values(['Order', 'L (SLP size)'])

    cols = [
        'Group (short)', 'Group (full)', 'Order', 'Family',
        'L (SLP size)', 'Seeds used',
        'CDS avg time/instr (ms)', 'CDS std (ms)', 'CDS CV (%)',
        'GAP avg time/instr (ms)', 'GAP std (ms)', 'GAP CV (%)',
        'Speedup (GAP/CDS)'
    ]

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        pd.DataFrame([['Per-L benchmark results']]).to_excel(
            writer, sheet_name='Plot2_PerLSize', index=False, header=False)
        df[cols].to_excel(
            writer, sheet_name='Plot2_PerLSize',
            startrow=2, index=False)

    n_groups = df['Group'].nunique()
    n_Ls     = df['L (SLP size)'].nunique()
    print(f"  saved  {Path(out_path).name}  ({n_groups} groups × {n_Ls} L values)")


# ── Plot 3: grand average, one row per group ──────────────────────────────────

def make_plot3(cpp, gap, out_path):
    # For each (group, L), average over seeds
    cpp_by_L = (cpp.groupby(['Group', 'Order', 'L'])['tpi']
                   .mean().reset_index())
    gap_by_L = (gap.groupby(['Group', 'L'])['tpi']
                   .mean().reset_index())

    # Then average the per-L means over all L values
    cpp_grand = (cpp_by_L.groupby(['Group', 'Order'])['tpi']
                         .agg(['mean', 'std', 'min', 'max'])
                         .reset_index())
    cpp_grand.columns = ['Group', 'Order',
                         'CDS grand avg time/instr (ms)',
                         'CDS std across L-means (ms)',
                         'CDS L-mean min (ms)', 'CDS L-mean max (ms)']
    cpp_grand['CDS L-sensitivity CV (%)'] = (
        cpp_grand['CDS std across L-means (ms)'] /
        cpp_grand['CDS grand avg time/instr (ms)'] * 100
    )

    gap_grand = (gap_by_L.groupby('Group')['tpi']
                         .agg(['mean', 'std', 'min', 'max'])
                         .reset_index())
    gap_grand.columns = ['Group',
                         'GAP grand avg time/instr (ms)',
                         'GAP std across L-means (ms)',
                         'GAP L-mean min (ms)', 'GAP L-mean max (ms)']
    gap_grand['GAP L-sensitivity CV (%)'] = (
        gap_grand['GAP std across L-means (ms)'] /
        gap_grand['GAP grand avg time/instr (ms)'] * 100
    )

    df = pd.merge(cpp_grand, gap_grand, on='Group')
    df['Group (short)'] = df['Group'].apply(short_name)
    df['Group (full)']  = df['Group']
    df['Family']        = df['Group'].apply(family)

    n_Ls   = cpp['L'].nunique()
    n_seed = cpp.groupby('Group')['Seed'].nunique().max()
    L_min  = int(cpp['L'].min())
    L_max  = int(cpp['L'].max())
    df['L sizes averaged'] = n_Ls
    df['Seeds per L']      = n_seed
    df['L range']          = f"{L_min}-{L_max}"
    df['Speedup (GAP/CDS)'] = (
        df['GAP grand avg time/instr (ms)'] /
        df['CDS grand avg time/instr (ms)']
    )
    df = df.sort_values('Order')

    cols = [
        'Group (short)', 'Group (full)', 'Order', 'Family',
        'L sizes averaged', 'Seeds per L', 'L range',
        'CDS grand avg time/instr (ms)', 'CDS std across L-means (ms)',
        'CDS L-sensitivity CV (%)', 'CDS L-mean min (ms)', 'CDS L-mean max (ms)',
        'GAP grand avg time/instr (ms)', 'GAP std across L-means (ms)',
        'GAP L-sensitivity CV (%)', 'GAP L-mean min (ms)', 'GAP L-mean max (ms)',
        'Speedup (GAP/CDS)'
    ]

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        pd.DataFrame([['Grand average benchmark results']]).to_excel(
            writer, sheet_name='Plot3_GrandAverage', index=False, header=False)
        df[cols].to_excel(
            writer, sheet_name='Plot3_GrandAverage',
            startrow=2, index=False)

    print(f"  saved  {Path(out_path).name}  ({len(df)} groups, "
          f"{n_Ls} L values × {n_seed} seeds)")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap", default="benchmark_results.csv",
                    help="GAP timing CSV")
    ap.add_argument("--cpp", default="benchmark_results_c++.csv",
                    help="CDS timing CSV")
    ap.add_argument("--out", default="benchmark_results",
                    help="Output base name (files will be <out>_plot1_FixedL.xlsx etc.)")
    args = ap.parse_args()

    if not Path(args.gap).exists():
        print(f"ERROR: {args.gap} not found"); return
    if not Path(args.cpp).exists():
        print(f"ERROR: {args.cpp} not found"); return

    print(f"Reading CSVs: {args.gap}, {args.cpp}")
    cpp, gap = load(args.gap, args.cpp)
    print(f"  CDS rows: {len(cpp)},  GAP rows: {len(gap)}")
    print(f"  Groups: {cpp['Group'].nunique()},  "
          f"L values: {cpp['L'].nunique()},  "
          f"Seeds per group: {cpp.groupby('Group')['Seed'].nunique().max()}\n")

    print("Generating Excel files...")
    make_plot1(cpp, gap, f"{args.out}_plot1_FixedL.xlsx")
    make_plot2(cpp, gap, f"{args.out}_plot2_PerLSize.xlsx")
    make_plot3(cpp, gap, f"{args.out}_plot3_GrandAverage.xlsx")

    print("\nDone. Run plot_6.py or plot_all.sh to generate plots.")


if __name__ == "__main__":
    main()
