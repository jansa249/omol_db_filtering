#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor

# --- CONFIGURATION ---
CSV_INPUT = Path('./output_filtered_data/molecule_index.csv')
OUTPUT_GRID = Path('multiplets_grid.png')
MAX_MOLS_TO_PLOT = 100

def analyze_smiles_fragments(csv_path):
    # 1. Load and Clean
    if not csv_path.exists():
        print(f"Error: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path, sep=';')
    df.columns = df.columns.str.strip()

    # Function to count components (dots + 1)
    def count_parts(smi):
        if pd.isna(smi) or str(smi) == "SMILES_ERROR":
            return 0
        return str(smi).count('.') + 1

    # 2. Apply counting
    # Note: using the column names 'smiles1' and 'smiles2' from your snippet
    df['parts_1'] = df['smiles1'].apply(count_parts)
    df['parts_2'] = df['smiles2'].apply(count_parts)
    df['total_parts'] = df['parts_1'] + df['parts_2']

    print("--- Fragment Analysis ---")
    print(f"Total entries analyzed: {len(df)}")
    print("\nBreakdown of total components per row:")
    print(df['total_parts'].value_counts().sort_index())

    # 3. Filter for multiplets (> 2 components)
    multiplets = df[df['total_parts'] > 2]
    
    if not multiplets.empty:
        print(f"\nFound {len(multiplets)} entries with multiple components.")
        
        # Save CSV for record keeping
        multiplets.to_csv('investigate_multiplets.csv', sep=';', index=False)
        print("Detailed list saved to 'investigate_multiplets.csv'")

        # 4. Generate Grid Plot
        print(f"Generating grid plot for up to {MAX_MOLS_TO_PLOT} multiplets...")
        
        # Take a sample if there are too many
        sample_size = min(len(multiplets), MAX_MOLS_TO_PLOT)
        plot_df = multiplets.sample(n=sample_size, random_state=42)
        
        mols_to_draw = []
        legends = []

        for _, row in plot_df.iterrows():
            # Combine the two SMILES strings back into one for visualization
            combined_smi = f"{row['smiles1']}.{row['smiles2']}"
            mol = Chem.MolFromSmiles(combined_smi)
            
            if mol:
                rdDepictor.Compute2DCoords(mol)
                mols_to_draw.append(mol)
                # Legend shows ID and the part distribution (e.g., 1+2 parts)
                legends.append(f"ID: {row['id']}\n({row['parts_1']}+{row['parts_2']} parts)")

        if mols_to_draw:
            img = Draw.MolsToGridImage(
                mols_to_draw,
                molsPerRow=5,
                subImgSize=(400, 400),
                legends=legends
            )
            img.save(str(OUTPUT_GRID))
            print(f"✅ Grid plot saved to {OUTPUT_GRID}")
    else:
        print("\nNo multi-component entries found.")

if __name__ == '__main__':
    analyze_smiles_fragments(CSV_INPUT)
