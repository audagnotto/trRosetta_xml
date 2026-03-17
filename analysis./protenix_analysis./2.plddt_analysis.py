#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2. pLDDT Analysis
=================
Extract and visualise the per-residue predicted Local Distance Difference
Test (pLDDT) confidence scores stored in the B-factor column of Protenix
CIF output files.

Usage
-----
    python 2.plddt_analysis.py <predictions_dir/> [output_dir/]

Arguments
---------
    predictions_dir  directory with Protenix output CIF (or PDB) files
    output_dir       where to write plots and tables (default: current dir)

Outputs
-------
    plddt_per_residue.png    – per-residue pLDDT profile (mean ± SD)
    plddt_distribution.png   – violin/box plot of mean pLDDT per model
    plddt_heatmap.png        – heatmap of per-residue pLDDT across all models
    plddt_results.txt        – table of per-model mean pLDDT values
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_structure, get_plddt, get_prediction_files


# AlphaFold2/3 pLDDT colour bands (same as AF2 publication)
# Bands: very high (90-100), confident (70-90), low (50-70), very low (<50)
BAND_COLORS = ['#0053D6', '#65CBF3', '#FFDB13', '#FF7D45']
BAND_LIMITS = [90, 70, 50, 0]
BAND_LABELS = ['Very high (90–100)', 'Confident (70–90)',
               'Low (50–70)', 'Very low (<50)']


def _plddt_color(val):
    if val >= 90:
        return '#0053D6'
    if val >= 70:
        return '#65CBF3'
    if val >= 50:
        return '#FFDB13'
    return '#FF7D45'


def main(pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    pred_files = get_prediction_files(pred_dir)
    if not pred_files:
        print("No prediction files (.pdb / .cif) found in:", pred_dir)
        sys.exit(1)

    print(f"Predictions: {len(pred_files)} files in {pred_dir}\n")

    all_plddt = []       # (n_models, n_residues)
    mean_plddt = []      # one mean per model
    model_names = []
    common_res = None

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        struct = load_structure(pred_file, name)
        res_nums, plddt = get_plddt(struct)
        if not plddt:
            print(f"  WARNING: no pLDDT values found in {name} – skipping")
            continue
        all_plddt.append(np.array(plddt))
        mean_plddt.append(np.mean(plddt))
        model_names.append(os.path.splitext(name)[0])
        if common_res is None:
            common_res = res_nums
        print(f"  {name:50s}  mean pLDDT = {np.mean(plddt):6.1f}")

    if not all_plddt:
        print("No pLDDT data found.")
        sys.exit(1)

    mean_plddt = np.array(mean_plddt)
    min_len = min(len(p) for p in all_plddt)
    arr = np.array([p[:min_len] for p in all_plddt])   # (models, residues)
    residues = common_res[:min_len] if common_res else list(range(1, min_len + 1))

    mean_per_res = arr.mean(axis=0)
    std_per_res  = arr.std(axis=0)

    # ── Plot 1: Per-residue pLDDT profile ─────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    for row in arr:
        ax.plot(residues, row, color='steelblue', alpha=0.15, linewidth=0.7)
    ax.plot(residues, mean_per_res, color='steelblue', linewidth=2,
            label='Mean pLDDT')
    ax.fill_between(residues, mean_per_res - std_per_res,
                    mean_per_res + std_per_res, alpha=0.2, color='steelblue',
                    label='±1 SD')
    # pLDDT confidence bands
    for ymin, col, lbl in zip([90, 70, 50], BAND_COLORS[:3], BAND_LABELS[:3]):
        ax.axhline(ymin, color=col, linestyle='--', linewidth=0.9, alpha=0.7)
    ax.set_xlim(residues[0], residues[-1])
    ax.set_ylim(0, 100)
    ax.set_xlabel('Residue number', fontsize=12)
    ax.set_ylabel('pLDDT', fontsize=12)
    ax.set_title('Per-residue pLDDT: Protenix Predictions', fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'plddt_per_residue.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"\nSaved: {out1}")

    # ── Plot 2: Per-model mean pLDDT violin / strip plot ──────────────────
    fig, ax = plt.subplots(figsize=(max(5, len(model_names) * 0.5 + 2), 5))
    vp = ax.violinplot([arr[i] for i in range(len(arr))],
                       positions=range(len(arr)), showmedians=True)
    for body in vp['bodies']:
        body.set_facecolor('steelblue')
        body.set_alpha(0.6)
    # overlay mean dots
    for i, (mn, col_val) in enumerate(zip(mean_plddt,
                                          [_plddt_color(v) for v in mean_plddt])):
        ax.scatter(i, mn, color=col_val, zorder=5, s=60,
                   edgecolors='black', linewidth=0.5)
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_ylabel('pLDDT', fontsize=12)
    ax.set_title('Mean pLDDT Distribution per Model', fontsize=13)
    # horizontal guide lines
    for ymin, col in zip([90, 70, 50], BAND_COLORS[:3]):
        ax.axhline(ymin, color=col, linestyle='--', linewidth=0.9, alpha=0.7)
    plt.tight_layout()
    out2 = os.path.join(output_dir, 'plddt_distribution.png')
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"Saved: {out2}")

    # ── Plot 3: Heatmap of per-residue pLDDT across all models ────────────
    fig, ax = plt.subplots(figsize=(12, max(3, len(model_names) * 0.35 + 1)))
    im = ax.imshow(arr, aspect='auto', cmap='RdYlBu', vmin=0, vmax=100,
                   origin='upper',
                   extent=[residues[0], residues[-1], len(arr) - 0.5, -0.5])
    fig.colorbar(im, ax=ax, label='pLDDT')
    ax.set_yticks(range(len(model_names)))
    ax.set_yticklabels(model_names, fontsize=7)
    ax.set_xlabel('Residue number', fontsize=12)
    ax.set_ylabel('Model', fontsize=12)
    ax.set_title('Per-residue pLDDT Heatmap', fontsize=13)
    plt.tight_layout()
    out3 = os.path.join(output_dir, 'plddt_heatmap.png')
    plt.savefig(out3, dpi=300)
    plt.close()
    print(f"Saved: {out3}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'plddt_results.txt')
    with open(out_txt, 'w') as f:
        f.write("model\tmean_plddt\tmin_plddt\tmax_plddt\n")
        for name, row in zip(model_names, arr):
            f.write(f"{name}\t{row.mean():.2f}\t{row.min():.2f}\t{row.max():.2f}\n")
        f.write(f"\noverall_mean\t{arr.mean():.2f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nSummary:")
    print(f"  Overall mean pLDDT : {arr.mean():.1f}")
    best = int(np.argmax(mean_plddt))
    print(f"  Best model (pLDDT) : {model_names[best]}  ({mean_plddt[best]:.1f})")
    # confidence band breakdown
    counts = {lbl: int(np.sum((mean_per_res >= lo) & (mean_per_res < hi)))
              for lbl, (lo, hi) in zip(BAND_LABELS,
                                       [(90, 101), (70, 90), (50, 70), (0, 50)])}
    print("  Mean per-residue pLDDT band counts:")
    for lbl, cnt in counts.items():
        print(f"    {lbl}: {cnt} residues")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    pred_dir = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) > 2 else "."
    main(pred_dir, out_dir)
