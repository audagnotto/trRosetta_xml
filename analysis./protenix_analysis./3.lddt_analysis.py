#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3. LDDT Analysis
================
Calculate the Local Distance Difference Test (LDDT) score between each
Protenix prediction and the experimental reference structure, and compare
the per-residue LDDT with the Protenix pLDDT confidence scores.

LDDT definition (Mariani et al. 2013, Bioinformatics 29:2722)
-------------------------------------------------------------
For each residue i, consider all pairs (i, j) with j within an inclusion
radius R0 (default 15 Å) of residue i in the reference structure and with
|i – j| > 0 (excluding immediate neighbours, |i – j| ≤ 1, in some formulations
the default excludes nothing).
For each threshold τ ∈ {0.5, 1, 2, 4} Å check whether:
    |d_pred(i,j) – d_ref(i,j)| < τ
LDDT for residue i = mean over all thresholds of the fraction of preserved
distances.
Global LDDT = mean over all residues.

Usage
-----
    python 3.lddt_analysis.py <experimental.pdb> <predictions_dir/> [output_dir/]

Outputs
-------
    lddt_per_residue.png      – mean ± SD per-residue LDDT profile
    lddt_vs_plddt.png         – scatter of per-residue mean LDDT vs mean pLDDT
    lddt_global.png           – bar chart of per-model global LDDT
    lddt_results.txt          – table of per-model global LDDT scores
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from utils import (load_structure, get_ca_coords, get_plddt,
                   get_residue_numbers, get_prediction_files)

# LDDT parameters
INCLUSION_RADIUS = 15.0   # Å  – only pairs within this distance in reference
THRESHOLDS = [0.5, 1.0, 2.0, 4.0]   # Å


def calc_lddt(ref_coords, pred_coords, r0=INCLUSION_RADIUS,
              thresholds=None):
    """Compute per-residue LDDT between two Cα coordinate arrays.

    Parameters
    ----------
    ref_coords  : (N, 3) ndarray – reference Cα positions
    pred_coords : (M, 3) ndarray – predicted Cα positions (M ≥ N)
    r0          : float – inclusion radius in Å
    thresholds  : list of float – distance-difference thresholds in Å

    Returns
    -------
    per_res_lddt : (N,) ndarray  – LDDT per residue (0–1)
    global_lddt  : float
    """
    if thresholds is None:
        thresholds = THRESHOLDS

    n = min(len(ref_coords), len(pred_coords))
    ref  = ref_coords[:n]
    pred = pred_coords[:n]

    # Pairwise Cα–Cα distances in reference and prediction
    diff_ref  = ref[:, None, :]  - ref[None, :, :]     # (N, N, 3)
    diff_pred = pred[:, None, :] - pred[None, :, :]
    d_ref  = np.sqrt((diff_ref  ** 2).sum(axis=-1))    # (N, N)
    d_pred = np.sqrt((diff_pred ** 2).sum(axis=-1))

    # Inclusion mask: pairs within r0 in the reference (excluding self)
    mask = (d_ref < r0) & ~np.eye(n, dtype=bool)

    # Absolute distance-difference
    dd = np.abs(d_pred - d_ref)                         # (N, N)

    per_res_lddt = np.zeros(n)
    for i in range(n):
        row_mask = mask[i]
        if row_mask.sum() == 0:
            per_res_lddt[i] = 0.0
            continue
        preserved = []
        for tau in thresholds:
            preserved.append((dd[i, row_mask] < tau).mean())
        per_res_lddt[i] = np.mean(preserved)

    return per_res_lddt, float(per_res_lddt.mean())


def main(exp_file, pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    exp_struct   = load_structure(exp_file, "experimental")
    ref_coords   = get_ca_coords(exp_struct)
    ref_res_nums = get_residue_numbers(exp_struct)

    pred_files = get_prediction_files(pred_dir)
    if not pred_files:
        print("No prediction files found in:", pred_dir)
        sys.exit(1)

    print(f"Reference : {exp_file}  ({len(ref_coords)} Cα atoms)")
    print(f"Predictions: {len(pred_files)} files\n")

    all_per_res = []
    global_lddts = []
    all_plddt = []
    model_names = []

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        pred_struct  = load_structure(pred_file, name)
        pred_coords  = get_ca_coords(pred_struct)

        per_res, glob = calc_lddt(ref_coords, pred_coords)
        all_per_res.append(per_res)
        global_lddts.append(glob)

        _, plddt_vals = get_plddt(pred_struct)
        all_plddt.append(np.array(plddt_vals) / 100.0)  # normalise to 0–1

        model_names.append(os.path.splitext(name)[0])
        print(f"  {name:50s}  LDDT = {glob:.4f}")

    global_lddts = np.array(global_lddts)
    best_idx     = int(np.argmax(global_lddts))
    min_len      = min(len(p) for p in all_per_res)
    arr          = np.array([p[:min_len] for p in all_per_res])
    mean_per_res = arr.mean(axis=0)
    std_per_res  = arr.std(axis=0)
    residues     = ref_res_nums[:min_len] if ref_res_nums else list(range(1, min_len + 1))

    # ── Plot 1: Global LDDT bar chart ─────────────────────────────────────
    fig_w = max(6, len(model_names) * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    colors = ['#4CAF50' if i == best_idx else '#2196F3'
              for i in range(len(global_lddts))]
    ax.bar(range(len(global_lddts)), global_lddts, color=colors, alpha=0.85,
           edgecolor='white')
    ax.axhline(global_lddts.mean(), color='gray', linestyle='--', linewidth=1.2,
               label=f'Mean = {global_lddts.mean():.4f}')
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Model', fontsize=12)
    ax.set_ylabel('Global LDDT', fontsize=12)
    ax.set_title('Global LDDT: Protenix Predictions vs Experimental', fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'lddt_global.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"\nSaved: {out1}")

    # ── Plot 2: Per-residue LDDT profile ──────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    for row in arr:
        ax.plot(residues, row[:min_len], color='#4CAF50', alpha=0.15,
                linewidth=0.7)
    ax.plot(residues, mean_per_res, color='#4CAF50', linewidth=2,
            label='Mean LDDT')
    ax.fill_between(residues, mean_per_res - std_per_res,
                    mean_per_res + std_per_res, alpha=0.2, color='#4CAF50',
                    label='±1 SD')
    ax.set_xlim(residues[0], residues[-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel('Residue number', fontsize=12)
    ax.set_ylabel('LDDT', fontsize=12)
    ax.set_title('Per-residue LDDT: Protenix vs Experimental', fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out2 = os.path.join(output_dir, 'lddt_per_residue.png')
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"Saved: {out2}")

    # ── Plot 3: LDDT vs pLDDT scatter ─────────────────────────────────────
    if all_plddt and all(len(p) >= min_len for p in all_plddt):
        plddt_arr  = np.array([p[:min_len] for p in all_plddt])
        lddt_flat  = arr.flatten()
        plddt_flat = plddt_arr.flatten()

        res_lr = stats.linregress(plddt_flat, lddt_flat)
        r2 = res_lr.rvalue ** 2

        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(plddt_flat, lddt_flat, alpha=0.1, s=8, color='steelblue',
                   rasterized=True)
        x_line = np.linspace(0, 1, 100)
        ax.plot(x_line, res_lr.intercept + res_lr.slope * x_line,
                'r-', linewidth=1.5, label=f'R² = {r2:.3f}')
        ax.plot([0, 1], [0, 1], 'k--', linewidth=0.8, alpha=0.4,
                label='Perfect calibration')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel('pLDDT (normalised)', fontsize=12)
        ax.set_ylabel('LDDT', fontsize=12)
        ax.set_title('LDDT vs pLDDT Calibration', fontsize=13)
        ax.legend(fontsize=10)
        plt.tight_layout()
        out3 = os.path.join(output_dir, 'lddt_vs_plddt.png')
        plt.savefig(out3, dpi=300)
        plt.close()
        print(f"Saved: {out3}")
        print(f"  LDDT vs pLDDT Pearson R² = {r2:.4f}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'lddt_results.txt')
    with open(out_txt, 'w') as f:
        f.write("model\tglobal_lddt\n")
        for name, score in zip(model_names, global_lddts):
            f.write(f"{name}\t{score:.6f}\n")
        f.write(f"\nmean_lddt\t{global_lddts.mean():.6f}\n")
        f.write(f"best_model\t{model_names[best_idx]}\t{global_lddts[best_idx]:.6f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nSummary:")
    print(f"  Mean global LDDT : {global_lddts.mean():.4f} ± {global_lddts.std():.4f}")
    print(f"  Best model       : {model_names[best_idx]}  ({global_lddts[best_idx]:.4f})")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    exp_file = sys.argv[1]
    pred_dir = sys.argv[2]
    out_dir  = sys.argv[3] if len(sys.argv) > 3 else "."
    main(exp_file, pred_dir, out_dir)
