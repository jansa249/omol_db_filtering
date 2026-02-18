#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

# --- IMPORT SHARED DICTIONARY ---
try:
    from queries import QUERY_DICT
except ImportError:
    # Manual fallback if queries.py is missing
    QUERY_DICT = {'phosphate': 'P~([O,OH])~([O,OH])~([O,OH])~[O,OH]'}

# --- CONFIGURATION ---
HITS_DIR = Path('.')
FILE_PATTERN = './analysis_results/matches_*.csv'
OUTPUT_GRID = Path('./grid_no_phosphates.png')
MAX_GRID_MOLS = 100
MOLS_PER_ROW = 10

def audit_phosphates():
    all_files = sorted(HITS_DIR.glob(FILE_PATTERN))
    if not all_files:
        print("No matches_*.csv files found.")
        return

    phos_query = Chem.MolFromSmarts(QUERY_DICT['phosphate'])
    
    # Data structures for tracking
    unique_all = set()
    missing_p_details = {}  # {(part, id): (smiles, sub_type)}
    summary_data = []

    print(f"{'Substructure':<15} | {'Total':<8} | {'With P':<8} | {'% With P':<10}")
    print("-" * 50)

    def has_p(smi):
        if pd.isna(smi) or str(smi) == "SMILES_ERROR": return False
        mol = Chem.MolFromSmiles(str(smi))
        return mol.HasSubstructMatch(phos_query) if mol else False

    # 1. Process Files and Print Table
    for file_path in all_files:
        sub_name = file_path.stem.replace('matches_', '')
        try:
            df = pd.read_csv(file_path, sep=';')
            df.columns = df.columns.str.strip()
            
            total_count = len(df)
            # Boolean mask for presence of phosphate
            p_mask = df['smiles1'].apply(has_p) | df['smiles2'].apply(has_p)
            p_count = p_mask.sum()
            p_perc = (p_count / total_count * 100) if total_count > 0 else 0
            
            print(f"{sub_name:<15} | {total_count:<8} | {p_count:<8} | {p_perc:>8.1f}%")

            # Update unique sets and "Missing" collection
            for _, row in df.iterrows():
                key = (row['dataset_part'], row['id'])
                unique_all.add(f"{key[0]}_{key[1]}")
                
                # If this specific molecule lacks phosphate, add to our "wanted" list
                if not (has_p(row['smiles1']) or has_p(row['smiles2'])):
                    if key not in missing_p_details:
                        combined_smi = f"{row['smiles1']}.{row['smiles2']}"
                        missing_p_details[key] = (combined_smi, sub_name)

        except Exception as e:
            print(f"Error processing {file_path.name}: {e}")

    print("-" * 50)
    print(f"Unique molecules across all files: {len(unique_all)}")
    print(f"Unique molecules MISSING phosphate: {len(missing_p_details)}")
    print("-" * 50)

    # 2. Print IDs of those missing phosphate
    if missing_p_details:
        print("\nList of IDs missing Phosphate:")
        print(f"{'Dataset Part':<15} | {'ID':<10} | {'Subtype'}")
        for (part, mol_id), (smi, sub) in missing_p_details.items():
            print(f"{str(part):<15} | {str(mol_id):<10} | {sub}")

        # 3. Draw Grid of the "Missing" ones
        mols_to_draw = []
        legends = []
        keys_to_plot = list(missing_p_details.keys())[:MAX_GRID_MOLS]

        for key in keys_to_plot:
            smi, sub = missing_p_details[key]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                rdDepictor.Compute2DCoords(mol)
                mols_to_draw.append(mol)
                legends.append(f"{key[0]}:{key[1]}\n({sub})")

        if mols_to_draw:
            img = Draw.MolsToGridImage(
                mols_to_draw,
                molsPerRow=MOLS_PER_ROW,
                subImgSize=(400, 400),
                legends=legends
            )
            img.save(str(OUTPUT_GRID))
            print(f"\n✅ Grid of missing phosphates saved to {OUTPUT_GRID}")

if __name__ == '__main__':
    audit_phosphates()
