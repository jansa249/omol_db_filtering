#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from queries import QUERY_DICT

# --- CONFIGURATION ---
CSV_INPUT = Path('./output_filtered_data/molecule_index.csv')
OUTPUT_DIR = Path('./analysis_results')
OUTPUT_DIR.mkdir(exist_ok=True)

N_SAMPLES_VIZ = 100
PHOS_QUERY = Chem.MolFromSmarts(QUERY_DICT['phosphate'])

def run_pipeline():
    # 1. Load Data Once
    print(f"--- Loading {CSV_INPUT.name} ---")
    df = pd.read_csv(CSV_INPUT, sep=';')
    df.columns = df.columns.str.strip()
    
    summary_stats = []
    
    # 2. Loop through all Substructures in QUERY_DICT
    for sub_name, smarts in QUERY_DICT.items():
        if sub_name == 'phosphate': continue # Skip the checker itself
        
        print(f"--- Processing: {sub_name} ---")
        query_mol = Chem.MolFromSmarts(smarts)
        
        # Define matching logic
        def get_mol_info(smi):
            if pd.isna(smi) or str(smi) == "SMILES_ERROR": return None
            return Chem.MolFromSmiles(str(smi))

        # Perform Substructure Search
        # We store the RDKit Mols temporarily to avoid re-parsing for Phosphate/Viz
        matches = []
        for idx, row in df.iterrows():
            m1 = get_mol_info(row['smiles1'])
            m2 = get_mol_info(row['smiles2'])
            
            has_sub = (m1.HasSubstructMatch(query_mol) if m1 else False) or \
                      (m2.HasSubstructMatch(query_mol) if m2 else False)
            
            if has_sub:
                has_phos = (m1.HasSubstructMatch(PHOS_QUERY) if m1 else False) or \
                           (m2.HasSubstructMatch(PHOS_QUERY) if m2 else False)
                matches.append({
                    **row.to_dict(),
                    'has_phosphate': has_phos,
                    'combined_mol': Chem.MolFromSmiles(f"{row['smiles1']}.{row['smiles2']}")
                })

        if not matches:
            print(f"No matches found for {sub_name}.")
            continue

        match_df = pd.DataFrame(matches)
        
        # 3. Save Substructure CSV
        match_df.drop(columns=['combined_mol']).to_csv(OUTPUT_DIR / f"matches_{sub_name}.csv", sep=';', index=False)

        # 4. Generate Summary Stats
        p_count = match_df['has_phosphate'].sum()
        p_perc = (p_count / len(match_df)) * 100
        summary_stats.append({'Substructure': sub_name, 'Count': len(match_df), 'With Phos': p_count, '%': p_perc})

        # 5. Visualization (Grid with Highlighting)
        sample_size = min(len(match_df), N_SAMPLES_VIZ)
        viz_sample = match_df.sample(n=sample_size)
        
        mols_to_draw = []
        hl_atoms = []
        legends = []

        for _, row in viz_sample.iterrows():
            mol = row['combined_mol']
            if mol:
                rdDepictor.Compute2DCoords(mol)
                mols_to_draw.append(mol)
                hl_atoms.append(list(mol.GetSubstructMatch(query_mol)))
                legends.append(f"{row['dataset_part']}:{row['id']} {'(+P)' if row['has_phosphate'] else '(no P)'}")

        img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=10, subImgSize=(400,400), 
                                   legends=legends, highlightAtomLists=hl_atoms)
        img.save(str(OUTPUT_DIR / f"grid_{sub_name}.png"))

    # 6. Final Summary Print
    print("\n" + "="*50)
    print(f"{'Substructure':<15} | {'Count':<8} | {'Phos %':<8}")
    print("-" * 35)
    for s in summary_stats:
        print(f"{s['Substructure']:<15} | {s['Count']:<8} | {s['%']:>7.1f}%")
    print("="*50)

if __name__ == '__main__':
    run_pipeline()
