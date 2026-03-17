# Protenix Analysis Suite

Python scripts for evaluating [Protenix](https://github.com/bytedance/protenix) structure
predictions against an experimental reference structure.

---

## Overview

| Script | Analysis | Key outputs |
|--------|----------|-------------|
| `1.rmsd_analysis.py` | Global & per-residue Cα RMSD | `rmsd_global.png`, `rmsd_per_residue.png`, `global_rmsd_results.txt` |
| `2.plddt_analysis.py` | pLDDT confidence scores | `plddt_per_residue.png`, `plddt_distribution.png`, `plddt_heatmap.png`, `plddt_results.txt` |
| `3.lddt_analysis.py` | Local Distance Difference Test (LDDT) + calibration vs pLDDT | `lddt_global.png`, `lddt_per_residue.png`, `lddt_vs_plddt.png`, `lddt_results.txt` |
| `4.contact_map_comparison.py` | Contact map precision / recall / F1 | `contact_map_experimental.png`, `contact_map_comparison.png`, `contact_precision_recall.png`, `contact_map_results.txt` |
| `5.secondary_structure_comparison.py` | Secondary structure composition (H/E/C) | `ss_composition.png`, `ss_per_residue.png`, `ss_agreement_heatmap.png`, `ss_results.txt` |
| `6.tm_score_analysis.py` | TM-score ranking | `tm_score_bar.png`, `tm_score_results.txt` |

A shared utility module `utils.py` provides structure loading, Cα extraction,
superimposition, pLDDT extraction, and other common functions used by all scripts.

---

## Requirements

```
biopython >= 1.79
numpy
matplotlib
scipy
```

Install via pip:

```bash
pip install biopython numpy matplotlib scipy
```

The secondary structure script (`5.secondary_structure_comparison.py`) optionally uses
the `mkdssp` / `dssp` binary for secondary structure assignment.  If it is not installed
the script falls back to assigning all residues as coil and will still run without error.
Install DSSP with:

```bash
# Debian/Ubuntu
sudo apt-get install dssp

# conda
conda install -c salilab dssp
```

---

## Input data

| Input | Description |
|-------|-------------|
| `experimental.pdb` | Reference experimental structure (PDB or mmCIF format) |
| `<predictions_dir>/` | Directory containing all Protenix output files (`.cif` or `.pdb`) |

Protenix typically names its output files `seed-{N}_sample-{M}_model.cif`.
The scripts accept any mix of `.pdb` / `.cif` / `.mmcif` files in the predictions
directory.  **pLDDT** confidence values are read from the B-factor column of the
prediction files (Protenix stores pLDDT there, just like AlphaFold2/3).

---

## Usage

Each script is invoked as:

```bash
python <script.py>  <experimental_structure>  <predictions_dir/>  [output_dir/]
```

`output_dir` defaults to the current directory if not given.

### Examples

```bash
# 1. RMSD analysis
python 1.rmsd_analysis.py ../6y75_experimental.pdb ../protenix_preds/ results/

# 2. pLDDT analysis (no experimental structure needed)
python 2.plddt_analysis.py ../protenix_preds/ results/

# 3. LDDT analysis
python 3.lddt_analysis.py ../6y75_experimental.pdb ../protenix_preds/ results/

# 4. Contact map comparison
python 4.contact_map_comparison.py ../6y75_experimental.pdb ../protenix_preds/ results/

# 5. Secondary structure comparison
python 5.secondary_structure_comparison.py ../6y75_experimental.pdb ../protenix_preds/ results/

# 6. TM-score analysis
python 6.tm_score_analysis.py ../6y75_experimental.pdb ../protenix_preds/ results/
```

---

## Metric descriptions

### RMSD (script 1)
Root Mean Square Deviation of Cα atoms after optimal superimposition.  Lower is
better.  Reported in Å.

### pLDDT (script 2)
Per-residue predicted Local Distance Difference Test confidence score (0–100) as
output by Protenix and stored in the B-factor column.  Higher is more confident.

| pLDDT range | Confidence colour |
|-------------|-------------------|
| 90 – 100 | Very high (dark blue) |
| 70 – 90 | Confident (light blue) |
| 50 – 70 | Low (yellow) |
| 0 – 50 | Very low (orange) |

### LDDT (script 3)
Local Distance Difference Test measures how well inter-residue distances are
preserved within a 15 Å inclusion sphere.  Four distance-difference thresholds
(0.5, 1, 2, 4 Å) are averaged.  Ranges 0–1; higher is better.  The calibration
plot (LDDT vs pLDDT) shows how well pLDDT predicts actual accuracy.

### Contact map (script 4)
A contact is defined as two Cα atoms ≤ 8 Å apart with sequence separation ≥ 6.
Precision, recall, and F1-score measure how accurately predicted contacts match
the experimental contact map.

### Secondary structure (script 5)
Uses DSSP to assign H (helix), E (strand), or C (coil) to each residue.
Comparison includes composition bar chart, per-residue assignment strip, and
residue-level agreement heatmap.

### TM-score (script 6)
Template Modelling score based on the Zhang & Skolnick (2004) formulation.
Computed by iterative Cα superposition.  Ranges 0–1:

| TM-score | Interpretation |
|----------|----------------|
| > 0.5 | Same fold |
| > 0.7 | Strongly similar structures |
| ~ 1.0 | Near-identical structures |

---

## Output files per run

All output files are written to the specified `output_dir/`.  Text tables use
tab-separated format and can be opened directly in Excel / LibreOffice Calc.
