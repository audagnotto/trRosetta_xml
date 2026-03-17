#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4. Contact Map Comparison
=========================
Build and compare Cα contact maps for the experimental structure and each
Protenix prediction.  A contact is defined as two residues whose Cα atoms
are within a given distance cutoff (default 8 Å) with a minimum sequence
separation (default |i – j| ≥ 6, i.e. long-range contacts).

Metrics computed
----------------
  Precision  = TP / (TP + FP)
  Recall     = TP / (TP + FN)
  F1-score   = 2 · precision · recall / (precision + recall)

Usage
-----
    python 4.contact_map_comparison.py <experimental.pdb> <predictions_dir/> [output_dir/]

Outputs
-------
    contact_map_experimental.png   – reference contact map
    contact_map_comparison.png     – overlay / predicted vs reference heatmap
    contact_precision_recall.png   – bar chart of precision & recall per model
    contact_map_results.txt        – table of precision, recall, F1 per model
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_structure, get_ca_coords, get_prediction_files


DISTANCE_CUTOFF = 8.0    # Å   – contact threshold
MIN_SEQ_SEP     = 6      # residues  – minimum sequence separation


def build_contact_map(coords, cutoff=DISTANCE_CUTOFF, min_sep=MIN_SEQ_SEP):
    """Return a boolean symmetric contact matrix."""
    n = len(coords)
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    cm = dist < cutoff
    # zero out diagonal and near-diagonal (sequence separation < min_sep)
    for k in range(min_sep):
        np.fill_diagonal(cm[k:, :], False)
        np.fill_diagonal(cm[:, k:], False)
    np.fill_diagonal(cm, False)
    return cm


def precision_recall_f1(ref_cm, pred_cm):
    tp = int((ref_cm & pred_cm).sum())
    fp = int((~ref_cm & pred_cm).sum())
    fn = int((ref_cm & ~pred_cm).sum())
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)
    return precision, recall, f1


def main(exp_file, pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    exp_struct = load_structure(exp_file, "experimental")
    ref_coords = get_ca_coords(exp_struct)
    ref_cm     = build_contact_map(ref_coords)

    pred_files = get_prediction_files(pred_dir)
    if not pred_files:
        print("No prediction files found in:", pred_dir)
        sys.exit(1)

    print(f"Reference : {exp_file}  ({len(ref_coords)} Cα atoms, "
          f"{ref_cm.sum()//2} contacts)")
    print(f"Cutoff    : {DISTANCE_CUTOFF} Å,  min sequence separation: {MIN_SEQ_SEP}\n")

    # ── Plot 1: Experimental contact map ──────────────────────────────────
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(ref_cm, cmap='Greys', origin='lower', interpolation='nearest')
    ax.set_xlabel('Residue', fontsize=11)
    ax.set_ylabel('Residue', fontsize=11)
    ax.set_title(f'Experimental Contact Map (< {DISTANCE_CUTOFF} Å)', fontsize=12)
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'contact_map_experimental.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"Saved: {out1}")

    precisions = []
    recalls    = []
    f1s        = []
    model_names = []
    pred_cms   = []

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        pred_struct = load_structure(pred_file, name)
        pred_coords = get_ca_coords(pred_struct)

        n = min(len(ref_coords), len(pred_coords))
        pred_cm = build_contact_map(pred_coords[:n])
        ref_cm_n = build_contact_map(ref_coords[:n])

        prec, rec, f1 = precision_recall_f1(ref_cm_n, pred_cm)
        precisions.append(prec)
        recalls.append(rec)
        f1s.append(f1)
        pred_cms.append(pred_cm)
        model_names.append(os.path.splitext(name)[0])
        print(f"  {name:50s}  P={prec:.3f}  R={rec:.3f}  F1={f1:.3f}")

    # ── Plot 2: Comparison heatmap (best model overlay) ───────────────────
    best_idx = int(np.argmax(f1s))
    best_pred_cm = pred_cms[best_idx]
    n_res = min(len(ref_coords), len(best_pred_cm))
    ref_trim = build_contact_map(ref_coords[:n_res])

    # Colour-coded overlay: TP=green, FP=red, FN=orange, TN=white
    overlay = np.zeros((n_res, n_res, 3), dtype=float)
    overlay[:] = [1, 1, 1]                                       # TN white
    overlay[ref_trim & best_pred_cm]   = [0.18, 0.60, 0.18]     # TP green
    overlay[~ref_trim & best_pred_cm]  = [0.85, 0.20, 0.20]     # FP red
    overlay[ref_trim & ~best_pred_cm]  = [0.95, 0.55, 0.10]     # FN orange

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].imshow(ref_trim, cmap='Greys', origin='lower',
                   interpolation='nearest')
    axes[0].set_title('Experimental', fontsize=11)
    axes[0].set_xlabel('Residue')
    axes[0].set_ylabel('Residue')

    axes[1].imshow(overlay, origin='lower', interpolation='nearest')
    axes[1].set_title(f'Best prediction: {model_names[best_idx]}\n'
                      f'P={precisions[best_idx]:.3f}  '
                      f'R={recalls[best_idx]:.3f}  '
                      f'F1={f1s[best_idx]:.3f}', fontsize=10)
    axes[1].set_xlabel('Residue')
    axes[1].set_ylabel('Residue')
    legend_handles = [
        mpatches.Patch(color='#2E9932', label='TP (both)'),
        mpatches.Patch(color='#D93333', label='FP (pred only)'),
        mpatches.Patch(color='#F28C1A', label='FN (exp only)'),
        mpatches.Patch(color='white',   label='TN (neither)',
                       edgecolor='gray'),
    ]
    axes[1].legend(handles=legend_handles, fontsize=8,
                   loc='upper right', framealpha=0.8)
    plt.tight_layout()
    out2 = os.path.join(output_dir, 'contact_map_comparison.png')
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"Saved: {out2}")

    # ── Plot 3: Precision & recall bar chart ──────────────────────────────
    x = np.arange(len(model_names))
    fig_w = max(6, len(model_names) * 0.5 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    w = 0.28
    ax.bar(x - w, precisions, w, label='Precision', color='#2196F3', alpha=0.85)
    ax.bar(x,     recalls,    w, label='Recall',    color='#4CAF50', alpha=0.85)
    ax.bar(x + w, f1s,        w, label='F1-score',  color='#FF9800', alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title(f'Contact Map Precision / Recall / F1  (cutoff {DISTANCE_CUTOFF} Å, '
                 f'sep ≥ {MIN_SEQ_SEP})', fontsize=12)
    ax.legend(fontsize=10)
    plt.tight_layout()
    out3 = os.path.join(output_dir, 'contact_precision_recall.png')
    plt.savefig(out3, dpi=300)
    plt.close()
    print(f"Saved: {out3}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'contact_map_results.txt')
    with open(out_txt, 'w') as fh:
        fh.write("model\tprecision\trecall\tf1_score\n")
        for name, p, r, f1 in zip(model_names, precisions, recalls, f1s):
            fh.write(f"{name}\t{p:.4f}\t{r:.4f}\t{f1:.4f}\n")
        fh.write(f"\nbest_model_f1\t{model_names[best_idx]}\t{f1s[best_idx]:.4f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nSummary (mean ± SD across {len(model_names)} models):")
    print(f"  Precision : {np.mean(precisions):.3f} ± {np.std(precisions):.3f}")
    print(f"  Recall    : {np.mean(recalls):.3f} ± {np.std(recalls):.3f}")
    print(f"  F1-score  : {np.mean(f1s):.3f} ± {np.std(f1s):.3f}")
    print(f"  Best model: {model_names[best_idx]}  (F1 = {f1s[best_idx]:.3f})")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    exp_file = sys.argv[1]
    pred_dir = sys.argv[2]
    out_dir  = sys.argv[3] if len(sys.argv) > 3 else "."
    main(exp_file, pred_dir, out_dir)
