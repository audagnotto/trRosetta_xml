#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
6. TM-score Analysis
====================
Calculate the Template Modelling score (TM-score) between each Protenix
prediction and the experimental reference structure.

TM-score (Zhang & Skolnick 2004, Proteins 57:702)
--------------------------------------------------
TM-score = (1 / L_target) * max_{rotation} Σ_i [ 1 / (1 + (d_i / d_0)²) ]

where:
  L_target  = number of residues in the reference structure
  d_i       = Cα distance between aligned residue pair i after superposition
  d_0(L)    = 1.24 * (L_target − 15)^(1/3) − 1.8   (d_0 ≥ 0.5 Å)

The maximisation over rotations is performed iteratively: the algorithm
starts from the superposition of all Cα pairs, then refines by keeping
only the subset of closest pairs, repeating until convergence.

TM-score interpretation
-----------------------
  > 0.5  proteins share the same fold
  > 0.7  strong similarity
  ~ 1.0  near-identical

Usage
-----
    python 6.tm_score_analysis.py <experimental.pdb> <predictions_dir/> [output_dir/]

Outputs
-------
    tm_score_bar.png        – bar chart of TM-scores
    tm_score_results.txt    – table of per-model TM-scores
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import Superimposer

sys.path.insert(0, os.path.dirname(__file__))
from utils import (load_structure, get_ca_atoms, get_ca_coords,
                   get_prediction_files)

# Fraction of closest-distance residues kept during each refinement iteration.
# 50 % is the standard TM-score iterative convergence setting.
CONVERGENCE_PERCENTILE = 50


def _d0(l_target):
    """Normalisation distance d0 as a function of chain length."""
    if l_target <= 21:
        return 0.5
    return max(0.5, 1.24 * (l_target - 15) ** (1.0 / 3.0) - 1.8)


def _tm_score_from_coords(ref, mob, d0_val):
    """Given two (N, 3) coordinate arrays, compute the TM-score contribution.

    Both arrays must already be optimally superimposed.
    """
    diff = ref - mob
    d2   = (diff ** 2).sum(axis=1)
    return (1.0 / (1.0 + d2 / d0_val**2)).sum()


def tm_score(ref_coords, mob_coords, max_iter=20):
    """Compute the TM-score by iterative superposition.

    Returns
    -------
    float  – TM-score in (0, 1]
    """
    n_ref = len(ref_coords)
    n_mob = len(mob_coords)
    n     = min(n_ref, n_mob)
    ref   = ref_coords[:n].copy()
    mob   = mob_coords[:n].copy()
    d0_val = _d0(n_ref)

    # Initial superposition using all Cα pairs
    def superpose(r, m):
        """In-place superposition of m onto r; returns aligned m."""
        r_center = r.mean(axis=0)
        m_center = m.mean(axis=0)
        r_c = r - r_center
        m_c = m - m_center
        H   = m_c.T @ r_c
        U, S, Vt = np.linalg.svd(H)
        det = np.linalg.det(Vt.T @ U.T)
        D   = np.diag([1, 1, det])
        R   = Vt.T @ D @ U.T
        m_aligned = (m_c @ R.T) + r_center
        return m_aligned, R, r_center, m_center

    mob_cur = mob.copy()
    best_tm  = 0.0
    best_mob = mob_cur.copy()

    for _ in range(max_iter):
        mob_aligned, R, r_ctr, m_ctr = superpose(ref, mob_cur)

        tm_val = _tm_score_from_coords(ref, mob_aligned, d0_val) / n_ref
        if tm_val > best_tm:
            best_tm  = tm_val
            best_mob = mob_aligned.copy()

        # Select residues with small distance for next iteration
        diff = ref - mob_aligned
        d    = np.sqrt((diff ** 2).sum(axis=1))
        d_thr = np.percentile(d, CONVERGENCE_PERCENTILE)
        mask  = d < d_thr
        if mask.sum() < 3:
            break

        ref_sub = ref[mask]
        mob_sub = mob_aligned[mask]
        mob_cur_sub, R2, _, _ = superpose(ref_sub, mob_sub)

        # Apply same rotation to all residues
        R_total = R2  # already from masked superpos
        m_ctr2  = mob_aligned[mask].mean(axis=0)
        r_ctr2  = ref[mask].mean(axis=0)
        mob_cur = ((mob_aligned - m_ctr2) @ R2.T) + r_ctr2

    return best_tm


def main(exp_file, pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    exp_struct = load_structure(exp_file, "experimental")
    ref_coords = get_ca_coords(exp_struct)
    l_target   = len(ref_coords)

    pred_files = get_prediction_files(pred_dir)
    if not pred_files:
        print("No prediction files found in:", pred_dir)
        sys.exit(1)

    print(f"Reference : {exp_file}  (L = {l_target},  d0 = {_d0(l_target):.3f} Å)")
    print(f"Predictions: {len(pred_files)} files\n")

    tm_scores  = []
    model_names = []

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        pred_struct = load_structure(pred_file, name)
        pred_coords = get_ca_coords(pred_struct)

        score = tm_score(ref_coords, pred_coords)
        tm_scores.append(score)
        model_names.append(os.path.splitext(name)[0])
        print(f"  {name:50s}  TM-score = {score:.4f}")

    tm_scores = np.array(tm_scores)
    best_idx  = int(np.argmax(tm_scores))

    # ── Plot: TM-score bar chart ───────────────────────────────────────────
    fig_w = max(6, len(model_names) * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    colors = ['#F44336' if i == best_idx else '#2196F3'
              for i in range(len(tm_scores))]
    ax.bar(range(len(tm_scores)), tm_scores, color=colors, alpha=0.85,
           edgecolor='white')
    ax.axhline(tm_scores.mean(), color='gray', linestyle='--', linewidth=1.2,
               label=f'Mean = {tm_scores.mean():.4f}')
    # Fold similarity thresholds
    ax.axhline(0.5, color='orange', linestyle=':', linewidth=1.0, alpha=0.8,
               label='Same fold (> 0.5)')
    ax.axhline(0.7, color='green',  linestyle=':', linewidth=1.0, alpha=0.8,
               label='Strong similarity (> 0.7)')
    ax.set_xticks(range(len(model_names)))
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Model', fontsize=12)
    ax.set_ylabel('TM-score', fontsize=12)
    ax.set_title('TM-score: Protenix Predictions vs Experimental', fontsize=13)
    ax.legend(fontsize=9)
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'tm_score_bar.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"\nSaved: {out1}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'tm_score_results.txt')
    with open(out_txt, 'w') as f:
        f.write(f"# d0 = {_d0(l_target):.4f} Å  (L_target = {l_target})\n")
        f.write("model\ttm_score\n")
        for name, score in zip(model_names, tm_scores):
            f.write(f"{name}\t{score:.6f}\n")
        f.write(f"\nmean_tm_score\t{tm_scores.mean():.6f}\n")
        f.write(f"best_model\t{model_names[best_idx]}\t{tm_scores[best_idx]:.6f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nSummary:")
    print(f"  Mean TM-score : {tm_scores.mean():.4f} ± {tm_scores.std():.4f}")
    print(f"  Best model    : {model_names[best_idx]}  ({tm_scores[best_idx]:.4f})")
    n_correct = int((tm_scores > 0.5).sum())
    print(f"  Models with TM-score > 0.5 (same fold): "
          f"{n_correct} / {len(tm_scores)}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    exp_file = sys.argv[1]
    pred_dir = sys.argv[2]
    out_dir  = sys.argv[3] if len(sys.argv) > 3 else "."
    main(exp_file, pred_dir, out_dir)
