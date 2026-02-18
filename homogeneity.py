#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from rdkit import Chem
from queries import QUERY_DICT

# --- CONFIGURATION ---
CSV_INPUT = Path('./output_filtered_data/molecule_index.csv')

# The names don't matter anymore, we just need the patterns
BASE_QUERIES = ['adenine', 'guanine', 'uracthym', 'cytosine']

def run_structural_bias_analysis():
    print(f"--- Loading {CSV_INPUT.name} ---")
    df = pd.read_csv(CSV_INPUT, sep=';')
    df.columns = df.columns.str.strip()

    # 1. Setup SMARTS (Your exact sugar strings)
    rna_smarts = 'C~C1~C(~O)~C(~O)~C~O1'
    dna_smarts = 'C~C1~C(~O)~C~C~O1'
    
    # Compile queries
    queries = {name: Chem.MolFromSmarts(smarts) for name, smarts in QUERY_DICT.items()}
    rna_q = Chem.MolFromSmarts(rna_smarts)
    dna_q = Chem.MolFromSmarts(dna_smarts)
    phos_q = queries.get('phosphate')
    
    # Group all base SMARTS into one list for a single check
    base_patterns = [queries[b] for b in BASE_QUERIES if b in queries]

    def classify_molecule(row):
        m1 = Chem.MolFromSmiles(str(row['smiles1']))
        m2 = Chem.MolFromSmiles(str(row['smiles2']))
        mols = [m for m in [m1, m2] if m]
        
        # Check for ANY base
        has_base = any(any(m.HasSubstructMatch(p) for p in base_patterns) for m in mols)
        
        # Check Sugar type
        is_rna = any(m.HasSubstructMatch(rna_q) for m in mols)
        is_dna = any(m.HasSubstructMatch(dna_q) for m in mols)
        
        # Check Phosphate
        has_p = any(m.HasSubstructMatch(phos_q) for m in mols)

        # High-level Classification
        sugar_label = "w/o sugar"
        if is_rna: sugar_label = "ribose"
        elif is_dna: sugar_label = "d-ribose"
            
        base_label = "w   base" if has_base else "w/o base"
        phos_label = "w   phos" if has_p else "w/o phos"

        profile = f"{sugar_label:>9} | {base_label:>8} | {phos_label:>8}"
        
        return pd.Series({
            'profile': profile,
            'is_full_nucleotide': (is_rna or is_dna) and has_base and has_p
        })

    print("--- Analyzing structural distributions... ---")
    analysis = df.apply(classify_molecule, axis=1)
    
    total = len(df)
    profile_counts = analysis['profile'].value_counts()
    nuc_count = analysis['is_full_nucleotide'].sum()

    print("\n" + "="*65)
    print("STRUCTURAL DIVERSITY REPORT")
    print("="*65)
    print(f"Total Unique Dimers: {total}")
    print("-" * 65)
    print(f"{'Structural Profile (Sugar | Base | Phos)':<45} | {'Count':<8} | {'%':<5}")
    print("-" * 65)
    
    for prof, count in profile_counts.items():
        perc = (count / total) * 100
        print(f"{prof:<45} | {count:<8} | {perc:>5.1f}%")

    print("-" * 65)
    print(f"TOTAL FULL NUCLEOTIDE BIAS: {(nuc_count/total)*100:>6.1f}%")
    print("="*65)

    if (nuc_count/total) > 0.95:
        print("\n[CONCLUSION] High structural redundancy.")
        print("The dataset consists almost entirely of complete nucleotides.")
        print("Model training will be biased toward complex, charged scaffolds.")

if __name__ == '__main__':
    run_structural_bias_analysis()
