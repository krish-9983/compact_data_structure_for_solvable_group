#!/usr/bin/env python3
"""
impact_plots.py — Generates impact2 and extra4 plots from raw CSVs.

Usage:
  python impact_plots.py <results_dir> <plots_dir>
"""

import sys, os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

results_dir = sys.argv[1]
plots_dir   = sys.argv[2]

GAP_COLS = [
    'Engine', 'Group', 'Order', 'Instructions',
    'Seed', 'Repetitions', 'Result', 'TotalTime_ms', 'TimePerRep_ms'
]

# ── Load CSVs ─────────────────────────────────────────────────────────────────

cpp = pd.read_csv(os.path.join(results_dir, "benchmark_results_c++.csv"))
cpp.columns = cpp.columns.str.strip()

gap = pd.read_csv(
    os.path.join(results_dir, "benchmark_results.csv"),
    header=None, names=GAP_COLS
)

info = pd.read_csv(os.path.join(results_dir, "group_info.csv"))
info.columns = info.columns.str.strip()

def fix_name(name):
    """Convert GAP group name format to standard format."""
    m = re.match(r'SmallGroup(\d+)_(\d+)$', name)
    if m:
        return f"SmallGroup({m.group(1)},{m.group(2)})"
    m = re.match(r'DirectProduct_?SmallGroup(\d+)_(\d+)_?CyclicGroup(\d+)', name)
    if m:
        return (f"DirectProduct(SmallGroup({m.group(1)},{m.group(2)}),"
                f"CyclicGroup({m.group(3)}))")
    return name

gap['Group'] = gap['Group'].apply(fix_name)

# ── Compute ns per operation ──────────────────────────────────────────────────

cpp['ns'] = cpp['TimePerRep_ms'] * 1e6 / cpp['Instructions']
gap['ns'] = gap['TimePerRep_ms'] * 1e6 / gap['Instructions']

# Grand average per group
cpp_avg = cpp.groupby('Group')['ns'].mean().reset_index().rename(columns={'ns': 'ns_cds'})
gap_avg = gap.groupby('Group')['ns'].mean().reset_index().rename(columns={'ns': 'ns_gap'})

# Start from group_info — it has all 100 groups with Order, Mode, BinSize, CayleySize
df = info[['Group', 'Order', 'Mode', 'BinSize', 'CayleySize']].copy()
df = pd.merge(df, cpp_avg, on='Group', how='left')
df = pd.merge(df, gap_avg, on='Group', how='left')
df = df.dropna(subset=['ns_cds', 'ns_gap'])
df = df.sort_values('Order')

print(f"Groups loaded: {len(df)}")
if len(df) != 100:
    print(f"  WARNING: expected 100 groups, got {len(df)}")

CDS = '#2196F3'
GAP = '#FF5722'
plt.style.use('seaborn-v0_8-whitegrid')

# ── Plot 4: grouped bar — avg time by mode ────────────────────────────────────

mode_avg = df.groupby('Mode')[['ns_cds', 'ns_gap']].mean()
modes    = mode_avg.index.tolist()
x = np.arange(len(modes))
w = 0.35

fig, ax = plt.subplots(figsize=(7, 5))
b1 = ax.bar(x - w/2, mode_avg['ns_cds'], w, label='CDS', color=CDS, edgecolor='white')
b2 = ax.bar(x + w/2, mode_avg['ns_gap'],  w, label='GAP', color=GAP, edgecolor='white')
ax.set_xticks(x)
ax.set_xticklabels([f'Mode {int(m)}' for m in modes], fontsize=12)
ax.set_ylabel('Avg per-operation time (ns)', fontsize=12)
ax.set_title('Average Query Time: CDS vs GAP by Construction Mode', fontsize=13)
ax.legend(fontsize=11)
for bar in b1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9)
for bar in b2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f'{bar.get_height():.0f}', ha='center', va='bottom', fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'impact2_grouped_bar_mode.png'), dpi=150)
plt.close()
print("  impact2_grouped_bar_mode.png done")

# ── Plot 5: space comparison — CDS vs Cayley table ───────────────────────────

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(df['Order'], df['CayleySize'].values / 1024,
        's-', color=GAP, label='Cayley table (n² × 4 bytes)', lw=2, ms=5)
ax.plot(df['Order'], df['BinSize'].values / 1024,
        'o-', color=CDS, label='CDS (precomputed.bin)', lw=2, ms=5)
ax.set_xlabel('Group Order', fontsize=12)
ax.set_ylabel('Space (KB)', fontsize=12)
ax.set_title('Space: CDS vs Cayley Table', fontsize=13)
ax.set_yscale('log')
ax.legend(fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'extra4_space_comparison.png'), dpi=150)
plt.close()
print("  extra4_space_comparison.png done")