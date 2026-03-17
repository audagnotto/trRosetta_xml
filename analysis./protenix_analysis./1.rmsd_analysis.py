#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1. RMSD Analysis
================
Calculate global and per-residue Cα RMSD between Protenix predictions
and an experimental reference structure.

Usage
-----
    python 1.rmsd_analysis.py <experimental.pdb> <predictions_dir/> [output_dir/]

Arguments
---------
    experimental   PDB or CIF file of the experimental reference structure
    predictions_dir  directory containing Protenix output CIF (or PDB) files
    output_dir       where to write plots and tables (default: current dir)

Outputs
-------
    rmsd_global.png           – bar chart of per-model global Cα RMSD
    rmsd_per_residue.png      – mean ± SD per-residue Cα RMSD profile
    global_rmsd_results.txt   – tab-separated table of model / RMSD values
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from utils import (load_structure, superimpose_ca, calc_per_residue_rmsd,
                   get_residue_numbers, get_prediction_files)


def main(exp_file, pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    exp_struct = load_structure(exp_file, "experimental")
    pred_files = get_prediction_files(pred_dir)

    if not pred_files:
        print("No prediction files (.pdb / .cif) found in:", pred_dir)
        sys.exit(1)

    print(f"Reference : {exp_file}")
    print(f"Predictions: {len(pred_files)} files in {pred_dir}\n")

    global_rmsds = []
    per_res_all = []
    model_names = []

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        pred_struct = load_structure(pred_file, name)

        # Global Cα RMSD (superimposition applied in-place)
        rmsd = superimpose_ca(exp_struct, pred_struct)
        global_rmsds.append(rmsd)

        # Per-residue Cα RMSD (re-loads mobile to avoid cumulative transforms)
        pred_struct2 = load_structure(pred_file, name + "_2")
        pr_rmsd = calc_per_residue_rmsd(exp_struct, pred_struct2)
        per_res_all.append(pr_rmsd)

        model_names.append(os.path.splitext(name)[0])
        print(f"  {name:50s}  Cα RMSD = {rmsd:6.2f} Å")

    global_rmsds = np.array(global_rmsds)
    best_idx = int(np.argmin(global_rmsds))

    # ── Plot 1: Global RMSD bar chart ──────────────────────────────────────
    fig_w = max(6, len(model_names) * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    colors = ['#2196F3' if i != best_idx else '#F44336'
              for i in range(len(global_rmsds))]
    ax.bar(range(len(global_rmsds)), global_rmsds, color=colors, alpha=0.85,
           edgecolor='white')
    ax.axhline(y=global_rmsds.mean(), color='gray', linestyle='--', linewidth=1.2,
               label=f'Mean = {global_rmsds.mean():.2f} Å')
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=8)
    ax.set_xlabel('Model', fontsize=12)
    ax.set_ylabel('Cα RMSD (Å)', fontsize=12)
    ax.set_title('Global Cα RMSD: Protenix Predictions vs Experimental', fontsize=13)
    ax.legend(fontsize=10)
    # annotate best
    ax.annotate(f'best\n{global_rmsds[best_idx]:.2f} Å',
                xy=(best_idx, global_rmsds[best_idx]),
                xytext=(best_idx + 0.5, global_rmsds[best_idx] + 0.1),
                fontsize=8, color='#F44336')
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'rmsd_global.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"\nSaved: {out1}")

    # ── Plot 2: Per-residue RMSD profile ──────────────────────────────────
    min_len = min(len(pr) for pr in per_res_all)
    arr = np.array([pr[:min_len] for pr in per_res_all])
    mean_pr = arr.mean(axis=0)
    std_pr  = arr.std(axis=0)

    ref_res = get_residue_numbers(exp_struct)[:min_len]
    residues = ref_res if ref_res else list(range(1, min_len + 1))

    fig, ax = plt.subplots(figsize=(10, 4))
    # individual traces
    for pr in arr:
        ax.plot(residues, pr[:min_len], color='steelblue', alpha=0.2,
                linewidth=0.7)
    ax.plot(residues, mean_pr, color='steelblue', linewidth=2,
            label='Mean RMSD')
    ax.fill_between(residues, mean_pr - std_pr, mean_pr + std_pr,
                    alpha=0.25, color='steelblue', label='±1 SD')
    ax.set_xlabel('Residue number', fontsize=12)
    ax.set_ylabel('Cα RMSD (Å)', fontsize=12)
    ax.set_title('Per-residue Cα RMSD: Protenix vs Experimental', fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out2 = os.path.join(output_dir, 'rmsd_per_residue.png')
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"Saved: {out2}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'global_rmsd_results.txt')
    with open(out_txt, 'w') as f:
        f.write("model\trmsd_angstrom\n")
        for name, rmsd in zip(model_names, global_rmsds):
            f.write(f"{name}\t{rmsd:.4f}\n")
        f.write(f"\nmean_rmsd\t{global_rmsds.mean():.4f}\n")
        f.write(f"std_rmsd\t{global_rmsds.std():.4f}\n")
        f.write(f"best_model\t{model_names[best_idx]}\t{global_rmsds[best_idx]:.4f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nSummary:")
    print(f"  Mean Cα RMSD : {global_rmsds.mean():.2f} ± {global_rmsds.std():.2f} Å")
    print(f"  Best model   : {model_names[best_idx]}  ({global_rmsds[best_idx]:.2f} Å)")
    print(f"  Worst model  : {model_names[int(np.argmax(global_rmsds))]}  ({global_rmsds.max():.2f} Å)")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    exp_file  = sys.argv[1]
    pred_dir  = sys.argv[2]
    out_dir   = sys.argv[3] if len(sys.argv) > 3 else "."
    main(exp_file, pred_dir, out_dir)
