#!/bin/bash
# ==============================================================================
# It generate the two impact plots used in the thesis:
#   impact2_grouped_bar_mode.png  — avg time CDS vs GAP by construction mode
#   extra4_space_comparison.png  — CDS vs Cayley table storage size
# ==============================================================================

set -uo pipefail

RESULTS_DIR="${1:-.}"
CPP_CSV="$RESULTS_DIR/benchmark_results_c++.csv"
GAP_CSV="$RESULTS_DIR/benchmark_results.csv"
INFO_CSV="$RESULTS_DIR/group_info.csv"
PLOTS_DIR="$RESULTS_DIR/plots"

echo "============================================================"
echo "PLOTTING — reading from $RESULTS_DIR"
echo "============================================================"

# Check required files exist
if [ ! -f "$CPP_CSV" ]; then echo "ERROR: $CPP_CSV not found"; exit 1; fi
if [ ! -f "$GAP_CSV" ]; then echo "ERROR: $GAP_CSV not found"; exit 1; fi
if [ ! -f "$INFO_CSV" ]; then echo "ERROR: $INFO_CSV not found"; exit 1; fi

mkdir -p "$PLOTS_DIR"

source $HOME/cds_env/bin/activate

# Run plot_6.py first to generate the three line plots
if [ -f plot_6.py ]; then
    echo "Running plot_6.py (line plots)..."
    cd "$RESULTS_DIR"
    python3 ../export_excel.py \
        --gap benchmark_results.csv \
        --cpp "benchmark_results_c++.csv" \
        --out benchmark_results
    python3 ../plot_6.py \
        --base benchmark_results \
        --dpi 150 \
        --out plots
    cd ..
    echo "  plot_6.py done"
fi

echo "Generating impact plots..."

python3 - "$RESULTS_DIR" "$PLOTS_DIR" << 'PYEOF'
import sys, os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

results_dir = sys.argv[1]
plots_dir   = sys.argv[2]

# ── Load data ─────────────────────────────────────────────────────────────────

GAP_COLS = [
    'Engine', 'Group', 'Order', 'Instructions',
    'Seed', 'Repetitions', 'Result', 'TotalTime_ms', 'TimePerRep_ms'
]

cpp = pd.read_csv(os.path.join(results_dir, "benchmark_results_c++.csv"))
cpp.columns = cpp.columns.str.strip()

gap = pd.read_csv(
    os.path.join(results_dir, "benchmark_results.csv"),
    header=None, names=GAP_COLS
)

def fix_gap_name(name):
    """Convert GAP group name format to standard format."""
    m = re.match(r'SmallGroup(\d+)_(\d+)$', name)
    if m:
        return f"SmallGroup({m.group(1)},{m.group(2)})"
    m = re.match(r'DirectProduct_?SmallGroup(\d+)_(\d+)_?CyclicGroup(\d+)', name)
    if m:
        return f"DirectProduct(SmallGroup({m.group(1)},{m.group(2)}),CyclicGroup({m.group(3)}))"
    return name

gap['Group'] = gap['Group'].apply(fix_gap_name)

info = pd.read_csv(os.path.join(results_dir, "group_info.csv"))
info.columns = info.columns.str.strip()

# Compute per-operation time in nanoseconds
cpp['ns'] = cpp['TimePerRep_ms'] * 1e6 / cpp['Instructions']
gap['ns'] = gap['TimePerRep_ms'] * 1e6 / gap['Instructions']

# Average per group
cpp_avg = cpp.groupby(['Group', 'Order', 'Mode'])['ns'].mean().reset_index()
gap_avg = gap.groupby(['Group', 'Order'])['ns'].mean().reset_index()

# Merge timing with group info
df = pd.merge(cpp_avg, gap_avg, on=['Group', 'Order'], suffixes=('_cds', '_gap'))
df = pd.merge(df, info[['Group', 'Mode', 'BinSize', 'CayleySize']], on='Group', how='left')
df['Mode'] = df['Mode_x'].fillna(df['Mode_y'])
df = df.sort_values('Order')

print(f"Groups loaded: {len(df)}")

# Colors
CDS = '#2196F3'
GAP = '#FF5722'
M3  = '#4CAF50'
M5  = '#FF9800'

plt.style.use('seaborn-v0_8-whitegrid')

# ── impact2: Grouped bar — avg time by construction mode ─────────────────────
mode_avg = df.groupby('Mode')[['ns_cds', 'ns_gap']].mean()
modes = mode_avg.index.tolist()
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
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            f'{bar.get_height():.1f}', ha='center', va='bottom', fontsize=9)
for bar in b2:
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
            f'{bar.get_height():.0f}', ha='center', va='bottom', fontsize=9)

fig.tight_layout()
fig.savefig(os.path.join(plots_dir, 'impact2_grouped_bar_mode.png'), dpi=150)
plt.close()
print("  impact2_grouped_bar_mode.png done")

# ── extra4: Space comparison — CDS vs Cayley table ───────────────────────────
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

PYEOF

echo ""
echo "============================================================"
echo "ALL PLOTS DONE — saved to $PLOTS_DIR/"
echo ""
echo "  From plot_6.py:"
echo "    plot1a_fixedL_vs_order.png"
echo "    plot2a_growingL_vs_order.png"
echo "    plot3a_grand_vs_order.png"
echo ""
echo "  Impact plots:"
echo "    impact2_grouped_bar_mode.png"
echo "    extra4_space_comparison.png"
echo "============================================================"