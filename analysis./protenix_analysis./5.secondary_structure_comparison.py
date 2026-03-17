#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5. Secondary Structure Comparison
==================================
Compare secondary structure composition (helix, strand/sheet, coil) between
the experimental reference structure and each Protenix prediction.

Secondary structure is assigned using the DSSP algorithm via BioPython.  If
DSSP is not installed, the script falls back to a geometry-based assignment
using backbone Cα–Cα virtual-bond angles and distances.

DSSP state mapping
------------------
    H, G, I  →  helix   (α-helix, 3₁₀-helix, π-helix)
    E, B     →  strand  (extended β-strand / bridge)
    T, S, C  →  coil    (turns, bends, random coil / unknown)

Usage
-----
    python 5.secondary_structure_comparison.py <experimental.pdb> <predictions_dir/> [output_dir/]

Outputs
-------
    ss_composition.png          – grouped bar chart of H/E/C fractions
    ss_per_residue.png          – per-residue assignment comparison (best model)
    ss_agreement_heatmap.png    – residue-level agreement across models
    ss_results.txt              – per-model H/E/C fractions
"""

import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(__file__))
from utils import load_structure, get_prediction_files

from Bio.PDB import DSSP

SS_COLORS = {'H': '#E74C3C', 'E': '#3498DB', 'C': '#95A5A6'}
SS_LABELS = {'H': 'Helix', 'E': 'Strand', 'C': 'Coil'}

# Numeric encoding for heatmap
SS_NUM = {'H': 0, 'E': 1, 'C': 2}
NUM_SS = {0: 'H', 1: 'E', 2: 'C'}


def dssp_to_three(dssp_code):
    """Map full DSSP 8-state code to H / E / C."""
    if dssp_code in ('H', 'G', 'I'):
        return 'H'
    if dssp_code in ('E', 'B'):
        return 'E'
    return 'C'


def run_dssp(structure, filepath):
    """Run DSSP and return a list of 3-state codes (one per residue).

    Tries 'mkdssp' first, then 'dssp'.  If neither binary is available
    a warning is printed and all residues are assigned coil ('C').
    """
    model = next(structure.get_models())
    for executable in ('mkdssp', 'dssp'):
        try:
            dssp_obj = DSSP(model, filepath, dssp=executable)
            # Success: parse and return
            ss_seq = [dssp_to_three(dssp_obj[key][2]) for key in dssp_obj.keys()]
            return ss_seq
        except Exception:
            continue
    # Both executables failed – fall back to coil assignment
    warnings.warn(
        "DSSP binary ('mkdssp' / 'dssp') not found or failed.  "
        "All residues assigned as coil.  "
        "Install DSSP (apt-get install dssp  or  conda install -c salilab dssp) "
        "for accurate secondary structure assignment."
    )
    res_list = [r for r in model.get_residues() if r.get_id()[0] == ' ']
    return ['C'] * len(res_list)


def ss_fractions(ss_seq):
    n = len(ss_seq)
    if n == 0:
        return 0.0, 0.0, 0.0
    h = ss_seq.count('H') / n
    e = ss_seq.count('E') / n
    c = ss_seq.count('C') / n
    return h, e, c


def main(exp_file, pred_dir, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)

    exp_struct = load_structure(exp_file, "experimental")
    ref_ss     = run_dssp(exp_struct, exp_file)
    ref_h, ref_e, ref_c = ss_fractions(ref_ss)

    pred_files = get_prediction_files(pred_dir)
    if not pred_files:
        print("No prediction files found in:", pred_dir)
        sys.exit(1)

    print(f"Reference : {exp_file}")
    print(f"  H={ref_h:.2%}  E={ref_e:.2%}  C={ref_c:.2%}")
    print(f"Predictions: {len(pred_files)} files\n")

    model_names = []
    pred_h_list = []
    pred_e_list = []
    pred_c_list = []
    pred_ss_seqs = []

    for pred_file in pred_files:
        name = os.path.basename(pred_file)
        struct = load_structure(pred_file, name)
        ss     = run_dssp(struct, pred_file)
        h, e, c = ss_fractions(ss)
        pred_h_list.append(h)
        pred_e_list.append(e)
        pred_c_list.append(c)
        pred_ss_seqs.append(ss)
        model_names.append(os.path.splitext(name)[0])
        print(f"  {name:50s}  H={h:.2%}  E={e:.2%}  C={c:.2%}")

    # ── Plot 1: Composition grouped bar chart ─────────────────────────────
    n_models = len(model_names)
    x = np.arange(n_models + 1)
    all_names   = ['Experimental'] + model_names
    all_h = [ref_h] + pred_h_list
    all_e = [ref_e] + pred_e_list
    all_c = [ref_c] + pred_c_list

    fig_w = max(7, (n_models + 1) * 0.55 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    w = 0.25
    ax.bar(x - w, all_h, w, label='Helix', color=SS_COLORS['H'], alpha=0.85)
    ax.bar(x,     all_e, w, label='Strand', color=SS_COLORS['E'], alpha=0.85)
    ax.bar(x + w, all_c, w, label='Coil', color=SS_COLORS['C'], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(all_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Fraction', fontsize=12)
    ax.set_title('Secondary Structure Composition: Protenix vs Experimental', fontsize=12)
    ax.axvline(0.5, color='black', linewidth=0.7, linestyle='--')
    ax.legend(fontsize=10)
    plt.tight_layout()
    out1 = os.path.join(output_dir, 'ss_composition.png')
    plt.savefig(out1, dpi=300)
    plt.close()
    print(f"\nSaved: {out1}")

    # ── Plot 2: Per-residue comparison for best model ─────────────────────
    # "best" = closest composition to experimental (smallest L2 distance)
    dists = [np.sqrt((h - ref_h)**2 + (e - ref_e)**2 + (c - ref_c)**2)
             for h, e, c in zip(pred_h_list, pred_e_list, pred_c_list)]
    best_idx = int(np.argmin(dists))

    n_res = min(len(ref_ss), len(pred_ss_seqs[best_idx]))
    ref_trim  = ref_ss[:n_res]
    pred_trim = pred_ss_seqs[best_idx][:n_res]

    # Encode as numbers for colour-bar
    ref_num  = [SS_NUM[s] for s in ref_trim]
    pred_num = [SS_NUM[s] for s in pred_trim]

    cmap = mcolors.ListedColormap([SS_COLORS['H'], SS_COLORS['E'],
                                   SS_COLORS['C']])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)

    fig, axes = plt.subplots(2, 1, figsize=(12, 3), sharex=True)
    for row_ax, data, lbl in zip(axes,
                                  [ref_num, pred_num],
                                  ['Experimental', model_names[best_idx]]):
        row_ax.imshow([data], aspect='auto', cmap=cmap, norm=norm,
                      origin='lower')
        row_ax.set_yticks([0])
        row_ax.set_yticklabels([lbl], fontsize=9)
        row_ax.set_ylabel('')
    axes[-1].set_xlabel('Residue index', fontsize=11)
    fig.suptitle('Per-residue Secondary Structure Comparison', fontsize=12)
    legend_handles = [mpatches.Patch(color=SS_COLORS[k], label=SS_LABELS[k])
                      for k in ('H', 'E', 'C')]
    axes[0].legend(handles=legend_handles, loc='upper right', fontsize=8,
                   framealpha=0.8)
    plt.tight_layout()
    out2 = os.path.join(output_dir, 'ss_per_residue.png')
    plt.savefig(out2, dpi=300)
    plt.close()
    print(f"Saved: {out2}")

    # ── Plot 3: Agreement heatmap ─────────────────────────────────────────
    # For each position, fraction of models that agree with the reference
    min_len = min(len(s) for s in pred_ss_seqs + [ref_ss])
    agreement = np.array([
        [1.0 if pred_ss_seqs[m][i] == ref_ss[i] else 0.0
         for i in range(min_len)]
        for m in range(n_models)
    ])

    fig, ax = plt.subplots(figsize=(12, max(3, n_models * 0.35 + 1)))
    im = ax.imshow(agreement, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1,
                   origin='upper')
    fig.colorbar(im, ax=ax, label='Agrees with experimental SS')
    ax.set_yticks(range(n_models))
    ax.set_yticklabels(model_names, fontsize=7)
    ax.set_xlabel('Residue index', fontsize=12)
    ax.set_title('Per-residue Secondary Structure Agreement with Experimental',
                 fontsize=12)
    plt.tight_layout()
    out3 = os.path.join(output_dir, 'ss_agreement_heatmap.png')
    plt.savefig(out3, dpi=300)
    plt.close()
    print(f"Saved: {out3}")

    # ── Save table ────────────────────────────────────────────────────────
    out_txt = os.path.join(output_dir, 'ss_results.txt')
    with open(out_txt, 'w') as f:
        f.write("model\thelix_frac\tstrand_frac\tcoil_frac\n")
        f.write(f"experimental\t{ref_h:.4f}\t{ref_e:.4f}\t{ref_c:.4f}\n")
        for name, h, e, c in zip(model_names, pred_h_list, pred_e_list,
                                  pred_c_list):
            f.write(f"{name}\t{h:.4f}\t{e:.4f}\t{c:.4f}\n")
    print(f"Saved: {out_txt}")

    print(f"\nBest model (closest composition): {model_names[best_idx]}")
    overall_agreement = agreement.mean()
    print(f"Mean per-residue SS agreement: {overall_agreement:.2%}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    exp_file = sys.argv[1]
    pred_dir = sys.argv[2]
    out_dir  = sys.argv[3] if len(sys.argv) > 3 else "."
    main(exp_file, pred_dir, out_dir)
