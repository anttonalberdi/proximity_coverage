#!/usr/bin/env python3
"""
Collect absolute paths to dereplicated genome files for a given assembly and binner.

Searches the standard coverage outputs:
1_single_coverage, 2_all_coverage, 3_proximity_coverage (per iteration),
and 4_random_coverage (per iteration) for the requested binner and assembly.
"""

import argparse
import csv
import glob
import os
from typing import Dict, Iterable, List


def gather_paths(assembly: str, binner: str) -> List[str]:
    """Return a de-duplicated, sorted list of genome file paths."""
    patterns = [
        f"1_single_coverage/{binner}_drep/{assembly}/dereplicated_genomes",
        f"2_all_coverage/{binner}_drep/{assembly}/dereplicated_genomes",
        f"3_proximity_coverage/*/{binner}_drep/{assembly}/dereplicated_genomes",
        f"4_random_coverage/*/{binner}_drep/{assembly}/dereplicated_genomes",
    ]

    files: List[str] = []
    seen = set()

    for pattern in patterns:
        for directory in glob.glob(pattern):
            if not os.path.isdir(directory):
                continue
            for path in sorted(glob.glob(os.path.join(directory, "*"))):
                if not os.path.isfile(path):
                    continue
                abspath = os.path.abspath(path)
                if abspath not in seen:
                    seen.add(abspath)
                    files.append(abspath)

    return files


def gather_genomeinfo_rows(assembly: str, binner: str) -> Dict[str, Dict[str, str]]:
    """Collect genomeInfo rows across runs keyed by genome filename."""
    patterns = [
        f"1_single_coverage/{binner}_drep/{assembly}/data_tables/genomeInfo.csv",
        f"2_all_coverage/{binner}_drep/{assembly}/data_tables/genomeInfo.csv",
        f"3_proximity_coverage/*/{binner}_drep/{assembly}/data_tables/genomeInfo.csv",
        f"4_random_coverage/*/{binner}_drep/{assembly}/data_tables/genomeInfo.csv",
    ]
    rows: Dict[str, Dict[str, str]] = {}
    for pattern in patterns:
        for path in glob.glob(pattern):
            if not os.path.isfile(path):
                continue
            with open(path, newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    genome = row.get("genome")
                    if genome:
                        rows.setdefault(genome, row)
    return rows


def write_paths(paths: Iterable[str], output: str) -> None:
    with open(output, "w") as handle:
        for p in paths:
            handle.write(f"{p}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="List dereplicated genome file paths for an assembly/binner."
    )
    parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        help="Assembly name (e.g., EHI01280)",
    )
    parser.add_argument(
        "-b",
        "--binner",
        required=True,
        choices=["metabat2", "maxbin2"],
        help="Binning workflow to query.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory to write outputs; created if missing (default: current directory).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    outfile = os.path.join(args.output_dir, f"{args.assembly}_{args.binner}_paths.txt")
    genomeinfo_out = os.path.join(args.output_dir, f"{args.assembly}_{args.binner}_genomeInfo.csv")
    paths = gather_paths(args.assembly, args.binner)
    write_paths(paths, outfile)

    info_rows = gather_genomeinfo_rows(args.assembly, args.binner)
    selected_genomes = {os.path.basename(p) for p in paths}
    filtered_rows = [info_rows[g] for g in selected_genomes if g in info_rows]

    if filtered_rows:
        fieldnames = ["genome", "completeness", "contamination", "length", "N50"]
        with open(genomeinfo_out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in filtered_rows:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
        print(f"Wrote {len(filtered_rows)} genomeInfo rows to {genomeinfo_out}")
    else:
        print("No matching genomeInfo entries found.")

    print(f"Wrote {len(paths)} paths to {outfile}")


if __name__ == "__main__":
    main()
