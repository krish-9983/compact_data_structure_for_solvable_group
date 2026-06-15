#!/usr/bin/env python3
"""
This generate CDS vs GAP line plots from the three Excel files.
Timings are in nanoseconds (ns).

Generates three plots used in the thesis:
  plot1a_fixedL_vs_order.png  — fixed L, time vs group order
  plot2a_growingL_vs_order.png — all L values, time vs group order
  plot3a_grand_vs_order.png   — grand average, time vs group order

Usage:
  python plot_6.py
  python plot_6.py --base benchmark_results --dpi 150 --out plots
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

warnings.filterwarnings("ignore")

matplotlib.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.color":        "#e0e0e0",
    "grid.linestyle":    "-",
    "grid.linewidth":    0.5,
    "axes.labelsize":    12,
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "legend.fontsize":   11,
    "legend.frameon":    True,
    "legend.framealpha": 0.9,
    "legend.edgecolor":  "none",
})

CDS_COL = "#1a8cd8"
GAP_COL = "#e24b4a"
MS_TO_NS = 1_000_000  # 1 ms = 1,000,000 ns


def read_excel(path):
    return pd.read_excel(path, skiprows=2)


def fmt_L(v):
    """Format L value: 1000 -> '1k'"""
    v = int(v)
    return f"{v // 1000}k" if v >= 1000 else str(v)


def save(fig, path, dpi):
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  {path.name}")


def draw_two_lines(ax, x, cds, gap, xlabel, ylabel, title,
                   log_x=False, xticks=None, xlabels=None):
    """Plot GAP (red) and CDS (blue) lines on ax."""
    ax.plot(x, gap, color=GAP_COL, lw=2.2, marker="s", ms=5, label="GAP")
    ax.plot(x, cds, color=CDS_COL, lw=2.2, marker="o", ms=5, label="CDS")
    if log_x:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: str(int(v)))
        )
    if xticks is not None:
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=45, ha="right")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()


# ── Plot 1a: fixed L, time vs group order ────────────────────────────────────
def plot_1a(f1, out_dir, dpi):
    """Time per operation at fixed L=1k, x = group order."""
    df = read_excel(f1)
    L    = int(df["L (fixed)"].iloc[0])
    nsds = int(df["Seeds used"].iloc[0])

    avg = (
        df.groupby("Order")
          .agg({
              "CDS avg time/instr (ms)": "mean",
              "GAP avg time/instr (ms)": "mean",
          })
          .reset_index()
          .sort_values("Order")
    )

    x   = avg["Order"].values
    cds = avg["CDS avg time/instr (ms)"].values * MS_TO_NS
    gap = avg["GAP avg time/instr (ms)"].values * MS_TO_NS

    fig, ax = plt.subplots(figsize=(9, 5))
    draw_two_lines(
        ax, x, cds, gap,
        xlabel="Group order  n",
        ylabel="Time per instruction  (ns)",
        title=f"Fixed L = {fmt_L(L)}  ·  avg over {nsds} random SLPs",
        log_x=True,
    )
    save(fig, out_dir / "plot1a_fixedL_vs_order.png", dpi)


# ── Plot 2a: all L values, time vs group order ────────────────────────────────
def plot_2a(f2, out_dir, dpi):
    """Time per operation for each L value, x = group order."""
    df   = read_excel(f2)
    Ls   = sorted(df["L (SLP size)"].unique())
    n_Ls = len(Ls)

    import matplotlib.cm as mcm
    cds_colors = [mcm.Blues(0.35 + 0.60 * i / max(n_Ls - 1, 1)) for i in range(n_Ls)]
    gap_colors = [mcm.Reds( 0.35 + 0.60 * i / max(n_Ls - 1, 1)) for i in range(n_Ls)]

    fig, ax = plt.subplots(figsize=(10, 5.5))

    for i, L in enumerate(Ls):
        sub = (
            df[df["L (SLP size)"] == L]
              .groupby("Order")
              .agg({
                  "CDS avg time/instr (ms)": "mean",
                  "GAP avg time/instr (ms)": "mean",
              })
              .reset_index()
              .sort_values("Order")
        )
        x   = sub["Order"].values
        cds = sub["CDS avg time/instr (ms)"].values * MS_TO_NS
        gap = sub["GAP avg time/instr (ms)"].values * MS_TO_NS
        lw  = 1.2 + 0.6 * i / max(n_Ls - 1, 1)
        lbl = fmt_L(L)
        ax.plot(x, cds, color=cds_colors[i], lw=lw, marker="o", ms=3,
                label=f"CDS L={lbl}")
        ax.plot(x, gap, color=gap_colors[i], lw=lw, marker="s", ms=3,
                ls="--", label=f"GAP L={lbl}", alpha=0.85)

    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: str(int(v))))
    ax.set_xlabel("Group order  n")
    ax.set_ylabel("Time per instruction  (ns)")
    ax.set_title(
        f"Growing L ({fmt_L(Ls[0])}→{fmt_L(Ls[-1])})  ·  each curve = avg over seeds\n"
        f"light→dark = small→large L  ·  x = group order"
    )
    from matplotlib.lines import Line2D
    handles = (
        [Line2D([0], [0], color=cds_colors[i], lw=1.5, marker="o", ms=4,
                label=f"CDS L={fmt_L(L)}") for i, L in enumerate(Ls)] +
        [Line2D([0], [0], color=gap_colors[i], lw=1.5, marker="s", ms=4,
                ls="--", label=f"GAP L={fmt_L(L)}") for i, L in enumerate(Ls)]
    )
    ax.legend(handles=handles, fontsize=7, ncol=4, framealpha=0.88)
    save(fig, out_dir / "plot2a_growingL_vs_order.png", dpi)


# ── Plot 3a: grand average, time vs group order ───────────────────────────────
def plot_3a(f3, out_dir, dpi):
    """Grand-average time per operation, x = group order."""
    df    = read_excel(f3)
    L_rng = str(df["L range"].iloc[0])
    n_Ls  = int(df["L sizes averaged"].iloc[0])
    n_s   = int(df["Seeds per L"].iloc[0])

    avg = (
        df.groupby("Order")
          .agg({
              "CDS grand avg time/instr (ms)": "mean",
              "GAP grand avg time/instr (ms)": "mean",
          })
          .reset_index()
          .sort_values("Order")
    )

    x   = avg["Order"].values
    cds = avg["CDS grand avg time/instr (ms)"].values * MS_TO_NS
    gap = avg["GAP grand avg time/instr (ms)"].values * MS_TO_NS

    fig, ax = plt.subplots(figsize=(9, 5))
    draw_two_lines(
        ax, x, cds, gap,
        xlabel="Group order  n",
        ylabel="Grand avg time per instruction  (ns)",
        title=f"Grand average  ·  avg over {n_Ls} L sizes ({L_rng}) × {n_s} seeds",
        log_x=True,
    )
    save(fig, out_dir / "plot3a_grand_vs_order.png", dpi)


# ── main ─────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="benchmark_results")
    ap.add_argument("--dpi",  type=int, default=150)
    ap.add_argument("--out",  default="plots")
    args = ap.parse_args()

    f1 = f"{args.base}_plot1_FixedL.xlsx"
    f2 = f"{args.base}_plot2_PerLSize.xlsx"
    f3 = f"{args.base}_plot3_GrandAverage.xlsx"

    for f in [f1, f2, f3]:
        if not Path(f).exists():
            print(f"ERROR: {f} not found — run export_excel.py first.")
            return

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output: {out_dir}/   Unit: nanoseconds (ns)\n")

    print("Plot 1a — fixed L, time vs group order:")
    plot_1a(f1, out_dir, args.dpi)

    print("Plot 2a — all L values, time vs group order:")
    plot_2a(f2, out_dir, args.dpi)

    print("Plot 3a — grand average, time vs group order:")
    plot_3a(f3, out_dir, args.dpi)

    print(f"\nDone. 3 plots saved in {out_dir}/")
    print("  plot1a_fixedL_vs_order.png")
    print("  plot2a_growingL_vs_order.png")
    print("  plot3a_grand_vs_order.png")


if __name__ == "__main__":
    main()