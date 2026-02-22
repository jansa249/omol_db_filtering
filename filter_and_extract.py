#!/usr/bin/env python3

import re
import csv
import tqdm
import numpy as np
from pathlib import Path
from scipy import sparse

from ase.db import connect
from ase.neighborlist import natural_cutoffs, NeighborList
from ase.io import write
from rdkit import Chem

# --- CONFIGURATION ---
CUTOFF_MULTIPLIER = 0.91                            # multiplier for 'natural_cutoffs', which determines bonding
INPUT_DIR = Path('./train_4M')                      # dir with .aselmdb files
OUTPUT_DIR = Path('./output_filtered_data')         # dir where coords and 'molecule_index.csv' will be written to
OUTPUT_XYZ_DIR = OUTPUT_DIR / 'coordinates'
CSV_LOG_PATH = OUTPUT_DIR / 'molecule_index.csv'


OUTPUT_XYZ_DIR.mkdir(exist_ok=True, parents=True)

# --- UTILS ---

def ase_to_smiles(atoms):
    """
    Constructs a SMILES string using ONLY heavy atoms and single bonds.
    """
    # 1. Initialize an editable RDKit molecule
    mol = Chem.RWMol()
    
    # 2. Add only Heavy Atoms (Atomic Number > 1)
    # Map old ASE index to new RDKit index
    ase_to_rd_idx = {}
    for i, atom in enumerate(atoms):
        if atom.number > 1:
            rd_idx = mol.AddAtom(Chem.Atom(int(atom.number)))
            ase_to_rd_idx[i] = rd_idx
            
    # 3. Use ASE NeighborList to find bonds
    cutoffs = [c * CUTOFF_MULTIPLIER for c in natural_cutoffs(atoms)]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=False)
    nl.update(atoms)
    
    # 4. Add Bonds to the RDKit molecule
    for i in ase_to_rd_idx.keys():
        indices, offsets = nl.get_neighbors(i)
        for j in indices:
            if j in ase_to_rd_idx:
                # Add a single bond if it doesn't already exist
                if not mol.GetBondBetweenAtoms(ase_to_rd_idx[i], ase_to_rd_idx[j]):
                    mol.AddBond(ase_to_rd_idx[i], ase_to_rd_idx[j], Chem.BondType.SINGLE)
    
    # 5. Finalize and convert to SMILES
    res_mol = mol.GetMol()
    
    # Sanitization for bad geometries
    try:
        res_mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(res_mol)
        return Chem.MolToSmiles(res_mol)
    except:
        return "SMILES_ERROR"

# --- MOLECULAR GRAPH CONSTRUCTION ---

def get_molecule_info(structure):
    """Returns component count, heavy atom counts, and atom indices for the ASE entry."""

    # 1. Get coordinates from entry
    atoms = structure.toatoms()

    # 2. Calculate NeighborList by distances and cutoffs
    cutoffs = [c * CUTOFF_MULTIPLIER for c in natural_cutoffs(atoms)]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    # 3. Connectivity matrix from NeighborList
    matrix = nl.get_connectivity_matrix()
    
    # 4. Calculate number of disconnected components and indices
    # which atom belongs to what part
    n_components, component_array = sparse.csgraph.connected_components(matrix)
    
    # 5. Filter for heavy atoms
    atomic_numbers = atoms.get_atomic_numbers()
    is_heavy = atomic_numbers > 1
    
    # 6. Generate output data
    mol_data = []
    for mol_id in range(n_components):
        indices = np.where(component_array == mol_id)[0]
        heavy_count = sum(is_heavy[indices])
        mol_data.append({
            'indices': indices.tolist(),
            'heavy_count': int(heavy_count),
            'atoms_obj': atoms[indices.tolist()]
        })
    
    return n_components, mol_data

# --- FILTERS ---

def filter_pipeline(structure, mol_info):
    n_components, mol_data = mol_info
    
    # 1. Filter by data_id
    if structure.data.get('data_id') != 'biomolecules':
        return False
    
    # 2. Filter by component size (e.g., dimer)
    if n_components != 2:
        return False
    
    # 3. Filter by heavy atoms per component
    if not all(m['heavy_count'] > 5 for m in mol_data):
        return False
        
    # 4. Filter by elements
    allowed = {'C', 'H', 'N', 'O', 'S', 'P'}
    # Composition string parsing
    present_elements = set(re.findall(r'[A-Z][a-z]?', structure.data.composition))
    if not present_elements.issubset(allowed):
        return False

    return True

# --- MAIN PROCESS ---

def main():
    input_paths = sorted(INPUT_DIR.glob('*.aselmdb'))
    
    # Open CSV for writing
    with open(CSV_LOG_PATH, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=';')
        csv_writer.writerow(['dataset_part', 'id', 'smiles1', 'smiles2'])

        # Go through every database part
        for input_path in input_paths:
            print(f'Processing {input_path.name}...')
            part_name = input_path.stem

            set_output_xyz_dir = OUTPUT_XYZ_DIR / part_name
            set_output_xyz_dir.mkdir(exist_ok=True)
            
            with connect(input_path) as db_in:
                for structure in tqdm.tqdm(db_in.select()):
                    # Get structural info once
                    mol_info = get_molecule_info(structure)
                    
                    if filter_pipeline(structure, mol_info):
                        n_components, mol_data = mol_info
                        atoms_complex = structure.toatoms()
                        
                        # 1. Create Directory for this specific ID
                        struct_dir = set_output_xyz_dir / f"{structure.id}_{structure.data.composition}"
                        struct_dir.mkdir(exist_ok=True)
                        
                        # 2. Save XYZ files (Complex and Components)
                        write(struct_dir / 'complex.xyz', atoms_complex)
                        
                        smiles_list = []
                        for i, m in enumerate(mol_data):
                            # Save individual xyz
                            write(struct_dir / f'component_{i}.xyz', m['atoms_obj'])
                            
                            # Generate SMILES
                            smi = ase_to_smiles(m['atoms_obj'])
                            smiles_list.append(smi)
                        
                        # 3. Write database part, id, and SMILES to CSV
                        csv_writer.writerow([part_name, structure.id, *smiles_list])

    print(f"Workflow Complete. Results in {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
