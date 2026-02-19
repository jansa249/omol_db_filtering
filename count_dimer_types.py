import os
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from queries import QUERY_DICT

# --- CONFIGURATION ---
TOP_ROWS = 7
NUM_SAMPLES = 4  # Number of random samples per category pair
OUTPUT_IMG_DIR = './analysis_results/category_galleries'
os.makedirs(OUTPUT_IMG_DIR, exist_ok=True)

# ... [Keep your classify_sugar and classify_smiles functions here] ...
# Convert strings to Mol objects once for efficiency
BASES = ['adenine', 'guanine', 'uracthym', 'cytosine']
BASE_MOLS = [Chem.MolFromSmarts(QUERY_DICT[base]) for base in BASES]

print(BASES, BASE_MOLS)

def classify_sugar(mol):
    if mol.HasSubstructMatch(Chem.MolFromSmarts(QUERY_DICT['ribose'])):
        return "r"
    if mol.HasSubstructMatch(Chem.MolFromSmarts(QUERY_DICT['dribose'])):
        return "d"
    return None

def classify_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "other"
    
    # "OR" Logic: Check if it matches ANY of the four base patterns
    has_base = any(mol.HasSubstructMatch(base_query) for base_query in BASE_MOLS)
    
    # Check for other components
    # (Using your existing logic for Sugar and Phosphate)
    sugar_type = classify_sugar(mol)
    has_phos = mol.HasSubstructMatch(Chem.MolFromSmarts(QUERY_DICT['phosphate']))

    if has_base and sugar_type and has_phos:
        return f'{sugar_type}NP'
    if has_base and sugar_type:
        return f'{sugar_type}N'
    if sugar_type and has_phos:
        return f'{sugar_type}P'

    if sugar_type:
        return f'{sugar_type}'
    if has_base:
        return f'B'
    if has_phos:
        return f'P'
    return "other"


# --- DATA PROCESSING ---
df = pd.read_csv('./output_filtered_data/molecule_index.csv', sep=';')

categories = ["other", "rNP", "dNP", "rN", "dN", "rP", "dP", "r", "d", "B", "P"]
matrix = pd.DataFrame(0, index=categories, columns=categories)

# Dictionary to store counts and samples symmetrically
# Key: tuple(sorted([type_a, type_b]))
# Value: {'count': int, 'samples': list}
dimer_store = {}
for idx, c1 in enumerate(categories):
    for c2 in categories[idx:]:
        dimer_store[tuple(sorted([c1, c2], key=lambda x: categories.index(x)))] = {'count': 0, 'samples': []}

print("Processing and sampling dimers...")
for _, row in df.iterrows():
    type_a = classify_smiles(row['smiles1'])
    type_b = classify_smiles(row['smiles2'])
    
    # 1. Update Symmetric Matrix
    matrix.loc[type_a, type_b] += 1
    if type_a != type_b:
        matrix.loc[type_b, type_a] += 1
    
    # 2. Reservoir Sampling Logic
    pair_key = tuple(sorted([type_a, type_b], key=lambda x: categories.index(x)))
    dimer_store[pair_key]['count'] += 1
    current_count = dimer_store[pair_key]['count']
    
    # If we have space, add it
    if len(dimer_store[pair_key]['samples']) < NUM_SAMPLES:
        dimer_store[pair_key]['samples'].append({
            'smiles': f"{row['smiles1']}.{row['smiles2']}",
            'part': row['dataset_part'],
            'id': row['id']
        })
    else:
        # If full, replace an existing sample with decreasing probability
        # This is a standard algorithm for random sampling from a stream
        r = random.randint(0, current_count - 1)
        if r < NUM_SAMPLES:
            dimer_store[pair_key]['samples'][r] = {
                'smiles': f"{row['smiles1']}.{row['smiles2']}",
                'part': row['dataset_part'],
                'id': row['id']
            }

# --- MATRIX VISUALIZATION ---
# [Your existing matrix plotting code here]
new_order = matrix["other"].sort_values(ascending=False).head(TOP_ROWS).index.tolist()
# ... (rest of the plot logic) ...

# --- DRAWING RANDOM SAMPLES ---
print("Generating random molecule galleries...")
for i_cat in new_order:
    for j_cat in new_order:
        pair_key = tuple(sorted([i_cat, j_cat], key=lambda x: categories.index(x)))
        sample_list = dimer_store[pair_key]['samples']
        
        if not sample_list: continue
            
        mols, labels = [], []
        for item in sample_list:
            mol = Chem.MolFromSmiles(item['smiles'])
            if mol:
                mols.append(mol)
                labels.append(f"{item['part']}:{item['id']}")

        if mols:
            img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(400, 400), legends=labels)
            img.save(f"{OUTPUT_IMG_DIR}/{i_cat}_vs_{j_cat}.png")

print("Done! Galleries represent a random distribution of each category.")
