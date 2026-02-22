import subprocess
import sys

# Workflow of analysis scripts
# Make sure you run 'filter_and_extract.py' first
WORKFLOW = [
    "analysis.py",
    "count_fragments.py",
    "distance_analysis.py",
    # "phosphate_summary.py",
    # "homogeneity.py",
    "count_dimer_types.py",
]

def run_analysis():
    for script in WORKFLOW:
        print(f"--- Running: {script} ---")
        try:
            # Executes each script using the current Python interpreter
            subprocess.run([sys.executable, script], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred in {script}. Halting pipeline.")
            sys.exit(1)
    
    print("\n✅ Full analysis workflow completed successfully.")

if __name__ == "__main__":
    run_analysis()
