# Biomolecular Structure Analysis Pipeline

## Overview
This project provides a complete workflow for processing, filtering, and analyzing large datasets of biomolecular structures (specifically nucleotide dimers) stored in ASE database files (`.aselmdb`).

The pipeline performs three main tasks:
1.  **Extraction & Conversion:** Reads raw 3D geometries, filters them based on composition/size, and converts them into 2D chemical representations (SMILES).
2.  **Structural Validation:** Analyzes connectivity and fragmentation to ensure data quality.
3.  **Chemical Analysis:** Classifies molecules by substructure (Nucleobases, Sugars, Phosphates) and generates statistical reports on dataset homogeneity.

**Note:** The raw database files (large binaries, `.xyz` directories) are excluded from this repository via `.gitignore` to maintain a lightweight codebase.

## Repository Structure

```text
.
├── filter_and_extract.py   # Core ETL script: ASE db -> Filter -> SMILES/XYZ
├── analysis.py             # Main analysis logic (substructure matching)
├── count_fragments.py      # Quality control: checks for fragmented molecules
├── homogeneity.py          # Statistical check for dataset bias
├── phosphate_summary.py    # Checks phosphate group presence
├── queries.py              # SMARTS patterns for chemical matching
├── environment.yml         # Conda environment configuration
└── output_filtered_data/   # (Generated) Stores processed results

## Setup

This project uses Conda for dependency management to ensure reproducibility with RDKit and ASE.

- Clone the repository:

```bash
git clone [https://github.com/yourusername/your-repo-name.git](https://github.com/yourusername/your-repo-name.git)
cd your-repo-name
```

- Create the environment:

```bash
conda env create -f environment.yml
```

- Activate the environment:

```bash
conda activate <env_name>
```

## Usage
### 1. Data Preparation

Place your raw .aselmdb database files into a directory named train_4M/ in the project root.
(Note: This directory is ignored by Git).
### 2. Extraction Pipeline

Run the extraction script to process the raw databases. This filters molecules (keeping only dimers with specific elements) and generates the index CSV.

```bash
python filter_and_extract.py
```
Input: train_4M/*.aselmdb
Output: output_filtered_data/molecule_index.csv and output_filtered_data/coordinates/

### 3. Analysis & Validation

Once the data is extracted, run the analysis tools in any order:
Substructure Search: Identifies nucleobases and generates grid visualizations.

```bash
python analysis.py
```

Homogeneity Check: Calculates the distribution of Sugar/Base/Phosphate profiles to detect dataset bias.

```bash
python homogeneity.py
```

Fragment Quality Control: Checks for molecules that have broken into unexpected fragments.
```bash
python count_fragments.py
```

## Outputs

The pipeline generates several key outputs in analysis_results/ and output_filtered_data/:
- molecule_index.csv: Master list of all processed IDs and their SMILES.
- matches_*.csv: Subsets of data matching specific queries (e.g., Adenine, Guanine).
- grid_*.png: Visual grids of molecular structures for quick inspection.
- complex.xyz / component_*.xyz: Extracted 3D coordinates for every filtered molecule.
