#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions shared across Protenix analysis scripts.

Provides helpers for loading PDB/CIF structures, extracting Cα atoms,
superimposing structures, and reading pLDDT values stored in B-factors.
"""

import os
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, Superimposer


def load_structure(filepath, structure_id="structure"):
    """Load a PDB or mmCIF file and return a BioPython Structure object.

    Parameters
    ----------
    filepath     : str – path to the structure file (.pdb, .cif, or .mmcif)
    structure_id : str – identifier assigned to the Structure object

    Returns
    -------
    Bio.PDB.Structure.Structure

    Raises
    ------
    FileNotFoundError if *filepath* does not exist.
    Bio.PDB.PDBExceptions.PDBConstructionException on malformed input.
    """
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ('.cif', '.mmcif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id, filepath)


def get_ca_atoms(structure, chain_id=None):
    """Return an ordered list of Cα atoms for standard amino-acid residues.

    Parameters
    ----------
    structure : Bio.PDB.Structure
    chain_id  : str or None – when None all chains are included;
                when a string only that chain is processed.

    Returns
    -------
    list of Bio.PDB.Atom  – one atom per standard residue that has a Cα atom
    """
    ca_atoms = []
    for model in structure:
        for chain in model:
            if chain_id and chain.get_id() != chain_id:
                continue
            for residue in chain:
                if PDB.is_aa(residue, standard=True) and 'CA' in residue:
                    ca_atoms.append(residue['CA'])
        break  # use only first MODEL
    return ca_atoms


def get_residue_numbers(structure, chain_id=None):
    """Return residue sequence numbers for standard amino-acid residues."""
    res_nums = []
    for model in structure:
        for chain in model:
            if chain_id and chain.get_id() != chain_id:
                continue
            for residue in chain:
                if PDB.is_aa(residue, standard=True) and 'CA' in residue:
                    res_nums.append(residue.get_id()[1])
        break
    return res_nums


def superimpose_ca(ref_structure, mob_structure, chain_id=None):
    """Superimpose *mob_structure* onto *ref_structure* using Cα atoms.

    **Modifies atom coordinates in mob_structure in place.**  Only the first
    MODEL of each structure is used.  Returns the post-superimposition RMSD (Å).
    """
    ref_ca = get_ca_atoms(ref_structure, chain_id)
    mob_ca = get_ca_atoms(mob_structure, chain_id)
    n = min(len(ref_ca), len(mob_ca))
    sup = Superimposer()
    sup.set_atoms(ref_ca[:n], mob_ca[:n])
    # Apply rotation/translation to all atoms in the mobile model
    for model in mob_structure:
        sup.apply(list(model.get_atoms()))
        break
    return sup.rms


def calc_per_residue_rmsd(ref_structure, mob_structure, chain_id=None):
    """Calculate per-residue Cα RMSD after global superposition.

    Returns
    -------
    list of float  – one value per aligned residue (Å)
    """
    ref_ca = get_ca_atoms(ref_structure, chain_id)
    mob_ca = get_ca_atoms(mob_structure, chain_id)
    n = min(len(ref_ca), len(mob_ca))
    ref_ca = ref_ca[:n]
    mob_ca = mob_ca[:n]

    sup = Superimposer()
    sup.set_atoms(ref_ca, mob_ca)
    for model in mob_structure:
        sup.apply(list(model.get_atoms()))
        break

    per_res = []
    for ra, ma in zip(ref_ca, mob_ca):
        diff = ra.get_vector() - ma.get_vector()
        per_res.append(diff.norm())
    return per_res


def get_plddt(structure, chain_id=None):
    """Extract pLDDT values from the B-factor column of a Protenix CIF.

    Returns
    -------
    res_nums : list of int
    plddt    : list of float  (0–100)
    """
    res_nums = []
    plddt = []
    for model in structure:
        for chain in model:
            if chain_id and chain.get_id() != chain_id:
                continue
            for residue in chain:
                if PDB.is_aa(residue, standard=True) and 'CA' in residue:
                    res_nums.append(residue.get_id()[1])
                    plddt.append(residue['CA'].get_bfactor())
        break
    return res_nums, plddt


def get_ca_coords(structure, chain_id=None):
    """Return a (N, 3) NumPy array of Cα coordinates."""
    ca_atoms = get_ca_atoms(structure, chain_id)
    return np.array([a.get_vector().get_array() for a in ca_atoms])


def get_prediction_files(pred_dir):
    """Return a sorted list of PDB/CIF files found in *pred_dir*.

    Parameters
    ----------
    pred_dir : str – path to the predictions directory

    Returns
    -------
    list of str – absolute paths to files with extensions
    .pdb / .cif / .mmcif (case-insensitive).  Returns an empty list if
    *pred_dir* does not exist or contains no supported files.
    """
    supported = ('.pdb', '.cif', '.mmcif')
    if not os.path.isdir(pred_dir):
        return []
    files = sorted([
        os.path.join(pred_dir, f)
        for f in os.listdir(pred_dir)
        if os.path.splitext(f)[1].lower() in supported
    ])
    return files
