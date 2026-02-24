import pandas as pd
from pathlib import Path
from ase.db import connect
from ase.io import write
from ase.neighborlist import neighbor_list
from collections import Counter
import tqdm
from datetime import datetime
import matplotlib.pyplot as plt

# --- Configuration ---
CSV_PATH = Path('./output_filtered_data/molecule_index.csv')
DATA_DIR = Path('./train_4M/')
OUTPUT_DIR = Path('./analysis_results/phosphate_structures')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH = Path("./analysis_results/phosphate_connectivity.txt")

def analyze_phosphate_environments(atoms):
    # Forced cutoffs for unoptimized P-O (2.1A), P-C (2.1A), P-H (1.6A)
    cutoffs = {
        ('P', 'O'): 2.1, ('P', 'C'): 2.1, ('P', 'H'): 1.6,
        ('O', 'H'): 1.3, ('O', 'C'): 1.7
    }
    
    i, j = neighbor_list('ij', atoms, cutoff=cutoffs)
    adj = {}
    for idx_i, idx_j in zip(i, j):
        adj.setdefault(idx_i, []).append(idx_j)
    
    p_indices = [a.index for a in atoms if a.symbol == 'P']
    p_data = []

    for p_idx in p_indices:
        p_neighbors = adj.get(p_idx, [])
        if not p_neighbors:
            p_data.append("Lone-P")
            continue

        groups = []
        for n_idx in p_neighbors:
            neighbor_atom = atoms[n_idx]
            if neighbor_atom.symbol == 'O':
                # Map what is on the Oxygen
                o_neighbors = sorted([atoms[idx].symbol for idx in adj.get(n_idx, []) if idx != p_idx])
                if not o_neighbors:
                    groups.append("(O)")
                else:
                    groups.append(f"(O-{'-'.join(o_neighbors)})")
            else:
                # Direct P-C, P-H, P-P
                groups.append(neighbor_atom.symbol)

        p_data.append("".join(sorted(groups)))
    return p_data

# --- Tracking ---
stats = {"looked_at": 0, "with_p": 0, "combos": []}
df = pd.read_csv(CSV_PATH, sep=';')

for part_name, group in tqdm.tqdm(df.groupby('dataset_part'), desc='Files', position=0):
    db_path = DATA_DIR / f"{part_name}.aselmdb"
    if not db_path.exists(): continue

    with connect(str(db_path)) as db:
        for _, row in tqdm.tqdm(group.iterrows(), total=len(group), desc=f'Structures in {part_name}', position=1, leave=False):
            stats["looked_at"] += 1
            if 'P' in str(row.get('smiles1', '')) or 'P' in str(row.get('smiles2', '')):
                try:
                    atoms = db.get(id=int(row['id'])).toatoms()
                    p_combos = analyze_phosphate_environments(atoms)
                    
                    if p_combos:
                        stats["with_p"] += 1
                        stats["combos"].extend(p_combos)
                        
                        # Multi-folder export
                        for p_type in set(p_combos):
                            target_dir = OUTPUT_DIR / p_type
                            target_dir.mkdir(exist_ok=True)
                            atoms.info['phosphate_types'] = "|".join(p_combos)
                            write(target_dir / f"{part_name}_{row['id']}.xyz", atoms, format='extxyz')
                except:
                    continue

# --- Generate and Save Report ---
summary_counts = Counter(stats["combos"])
report_content = [
    f"PHOSPHATE ANALYSIS REPORT - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    "="*60,
    f"Total structures scanned from CSV:  {stats['looked_at']}",
    f"Molecules containing P:            {stats['with_p']}",
    f"Total P centers analyzed:           {len(stats['combos'])}",
    "-"*60,
    f"{'COUNT':<10} | CONNECTIVITY TYPE",
    "-"*60
]

for combo, count in summary_counts.most_common():
    report_content.append(f"{count:<10} | {combo}")
report_content.append("="*60)

final_report = "\n".join(report_content)

# Print to console
print(f"\n{final_report}")

# Save to file
with open(SUMMARY_PATH, "w") as f:
    f.write(final_report)

print(f"Summary written to: {SUMMARY_PATH}")

# --- Visualization with Matplotlib ---
if stats["combos"]:
    # Prepare data (sorted by count from most_common)
    sorted_combos = summary_counts.most_common()
    types = [item[0] for item in sorted_combos][:10][::-1]
    counts = [item[1] for item in sorted_combos]
    total = sum(counts)
    counts = counts[:10][::-1]

    plt.figure(figsize=(6, len(types) * 0.4 + 2))
    
    # Using plain plt with some basic styling
    bars = plt.barh(types, counts, color=(1,0.49,0), edgecolor='white', height=0.8)
    
    plt.xscale('log')
    plt.xlabel('Occurrence Count (Log Scale)', fontsize=12, fontweight='bold')
    plt.ylabel('Connectivity Type', fontsize=12, fontweight='bold')
    plt.title(f'Distribution of P Environments\nN={total}', fontsize=14, pad=20)
    
    # Add light grid for log scale readability
    plt.grid(visible=True, which='both', axis='x', linestyle='--', alpha=0.5)
    plt.xlim(right=max(counts) * 8)

    # Add count labels next to bars
    for bar in bars:
        width = bar.get_width()
        plt.text(width * 1.1, bar.get_y() + bar.get_height()/2, 
                 f'{int(width):,}', va='center', fontsize=12)

    plt.tight_layout()

    # Save and Show
    plot_path = Path("./analysis_results/phosphorus_distribution.png")
    plt.savefig(plot_path, dpi=300)
    # plt.show()
    
    print(f"Distribution plot saved to: {plot_path}")
else:
    print("No Phosphorus types found to graph.")
