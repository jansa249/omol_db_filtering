#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shutil
from pathlib import Path
from ase.io import read
from scipy.spatial.distance import cdist
from tqdm import tqdm

# --- CONFIGURATION ---
XYZ_ROOT = Path('./output_filtered_data/coordinates')
OUTPUT_DIR = Path('./analysis_results')
CACHE_FILE = OUTPUT_DIR / 'distance_data_cache.csv'
C_OUTPUT_DIR = OUTPUT_DIR / 'distance_extremes'
OUTPUT_DIR.mkdir(exist_ok=True)

def get_min_dist(path_0, path_1):
    """Calculates min distance between heavy atoms of two components."""
    atoms0 = read(path_0)
    atoms1 = read(path_1)
    pos0 = atoms0.get_positions()[atoms0.get_atomic_numbers() > 1]
    pos1 = atoms1.get_positions()[atoms1.get_atomic_numbers() > 1]
    if len(pos0) == 0 or len(pos1) == 0:
        return np.nan
    return np.min(cdist(pos0, pos1))

def run_distance_analysis():
    # 1. LOAD CACHE OR RUN EXTRACTION
    if CACHE_FILE.exists():
        print(f"--- Loading cached data from {CACHE_FILE.name} ---")
        df = pd.read_csv(CACHE_FILE, sep=';')
    else:
        print("--- No cache found. Scanning XYZ Coordinates ---")
        results = []
        xyz_folders = list(XYZ_ROOT.glob('*/*'))
        for folder in tqdm(xyz_folders):
            comp0, comp1 = folder / 'component_0.xyz', folder / 'component_1.xyz'
            if comp0.exists() and comp1.exists():
                try:
                    part_name = folder.parent.name
                    struct_id = folder.name.split('_')[0]
                    results.append({
                        'dataset_part': part_name,
                        'id': struct_id,
                        'min_dist': get_min_dist(comp0, comp1),
                        'folder_path': str(folder)
                    })
                except Exception as e:
                    print(f"Error processing {folder}: {e}")
        df = pd.DataFrame(results).dropna(subset=['min_dist'])
        df.to_csv(CACHE_FILE, index=False, sep=';')

    df = df.sort_values('min_dist').reset_index(drop=True)

    # 2. GENERATE HISTOGRAM
    plt.figure(figsize=(5, 3))
    BINS = 50
    hist, bins = np.histogram(df['min_dist'], bins=BINS)

    # plt.bar((bins[1:]+bins[:-1])/2, np.log10(hist), width=np.max(df['min_dist'])/BINS, color=(1,0.49,0), edgecolor='white', alpha=0.8)
    # plt.ylabel('$\log_{10}(\mathrm{count})$')
    
    plt.bar((bins[1:]+bins[:-1])/2, hist, width=np.max(df['min_dist'])/BINS, color=(1,0.49,0), edgecolor='white', alpha=0.8)
    plt.ylabel('$\mathrm{count}$')

    plt.axvline(df['min_dist'].median(), color=(0.8,0,0), linestyle='--', label=f'Median: {df["min_dist"].median():.2f} Å')
    plt.xlabel('Distance [Å]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'distance_distribution.png', dpi=600)

    # 3. PRINT STATS
    print("\n" + "="*65)
    print(f"{'10 CLOSEST (POTENTIAL CLASHES)':^65}")
    print("-" * 65)
    print(df.head(10)[['dataset_part', 'id', 'min_dist']].to_string(index=False))
    print("\n" + "="*65)
    print(f"{'10 FARTHEST (LOOSE DIMERS)':^65}")
    print("-" * 65)
    print(df.tail(10)[['dataset_part', 'id', 'min_dist']].to_string(index=False))

    # 4. EXPORT EXTREME COORDINATES
    if C_OUTPUT_DIR.exists(): shutil.rmtree(C_OUTPUT_DIR)
    C_OUTPUT_DIR.mkdir(exist_ok=True)
    
    extremes = pd.concat([df.head(10), df.tail(10)])
    print(f"\n--- Exporting {len(extremes)} structures ---")
    
    for _, row in extremes.iterrows():
        label = "CLOSEST" if row['min_dist'] < df['min_dist'].median() else "FARTHEST"
        src_file = Path(row['folder_path']) / 'complex.xyz'
        if src_file.exists():
            new_name = f"{label}_{row['min_dist']:.2f}A_{row['dataset_part']}_{row['id']}.xyz"
            shutil.copy2(src_file, C_OUTPUT_DIR / new_name)
    
    print(f"✅ Coordinate files ready in: {C_OUTPUT_DIR}")

if __name__ == '__main__':
    run_distance_analysis()
