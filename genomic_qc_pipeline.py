"""Genomic Quality Control & Visualization Pipeline
Author: Irfan Chowdhury
Description: Processes, cleans, and merges call-rate and DNA concentration datasets, and generates QC box plots.
"""

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt

def run_qc_pipeline(call_rate_path: str, conc_path: str, output_plot: str, output_stats: str):
    """Processes QC datasets, merges them, generates a box plot, and exports summary statistics."""
    try:
        # Load call rates dataset
        call_rates = pd.read_csv(call_rate_path)
        call_rates = call_rates.set_index("sample_id")
        assert call_rates.index.is_unique, "Duplicate sample IDs found in call rate file."

        # Filter out samples with a call rate of 0 or less
        call_rates = call_rates.loc[call_rates["call_rate"] > 0, :]

        # Load concentration dataset (tab-separated)
        conc = pd.read_csv(conc_path, sep="\t")
        conc = conc.set_index("Sample ID")
        assert conc.index.is_unique, "Duplicate sample IDs found in concentration file."

        # Merge datasets on sample ID (index) and drop missing values
        data = call_rates.merge(conc, left_index=True, right_index=True)
        data = data.dropna()

        if data.empty:
            print("Error: The merged dataset is empty after filtering and dropping missing values.")
            sys.exit(1)

        # Find unique concentrations and sort them in increasing order
        unique_conc = data["DNA Concentration"].unique()
        unique_conc.sort()

        # Build vector of call rates per concentration
        call_rate_per_conc = [
            data.loc[data["DNA Concentration"] == c, "call_rate"]
            for c in unique_conc
        ]

        # Generate box plot
        plt.figure(figsize=(8, 6))
        n_samples = len(data)  
        plt.boxplot(x=call_rate_per_conc, tick_labels=unique_conc)
        plt.title(f"Call Rates vs DNA Concentrations (n={n_samples})", weight="bold")
        plt.xlabel("DNA Concentration (ng/µL)", weight="bold")
        plt.ylabel("Call Rate", weight="bold")
        plt.xticks(rotation=-20)
        plt.tight_layout()

        # Save figure
        plt.savefig(output_plot)
        print(f"QC box plot successfully saved to '{output_plot}'.")

        # Export summary statistics
        stats = data.groupby("DNA Concentration")["call_rate"].describe()
        stats.to_csv(output_stats)
        print(f"Summary statistics successfully saved to '{output_stats}'.")

    except FileNotFoundError as e:
        print(f"Error: A required file was not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error during QC pipeline execution: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Run genomic quality control pipeline and generate QC box plots.")
    parser.add_argument("-c", "--callrates", required=True, help="Path to call_rate.csv input file")
    parser.add_argument("-d", "--concentration", required=True, help="Path to concentration.tsv input file")
    parser.add_argument("-p", "--plot", default="boxplot.png", help="Filename for output boxplot (default: boxplot.png)")
    parser.add_argument("-s", "--stats", default="statistics.csv", help="Filename for output statistics CSV (default: statistics.csv)")

    args = parser.parse_args()
    run_qc_pipeline(args.callrates, args.concentration, args.plot, args.stats)

if __name__ == "__main__":
    main()