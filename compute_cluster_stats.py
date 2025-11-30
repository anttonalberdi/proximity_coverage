#!/usr/bin/env python3
"""
Compute clustering statistics from dRep Ndb and genomeInfo outputs.

Current functionality: report how many clusters contain at least one genome
from each binning type (e.g., MET_ALL, MET_R9).
"""

import argparse
import csv
import os
import re
from collections import defaultdict
from statistics import mean
from typing import Dict, Set


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute cluster statistics across binning types."
    )
    parser.add_argument(
        "-m",
        "--max-genomes-per-cluster",
        type=int,
        required=True,
        help="Maximum genomes expected per cluster (reserved for downstream stats).",
    )
    parser.add_argument(
        "-n",
        "--ndb",
        required=True,
        help="Path to Ndb.csv containing pairwise ANI and primary_cluster columns.",
    )
    parser.add_argument(
        "-g",
        "--genome-info",
        required=True,
        help="Path to genomeInfo.csv (currently parsed for future metrics).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=".",
        help="Directory to write outputs (created if missing).",
    )
    return parser.parse_args()


def binning_type(genome_name: str) -> str:
    """Extract binning type prefix (first two underscore-separated parts)."""
    base = os.path.basename(genome_name)
    parts = base.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return parts[0]


def load_cluster_members(ndb_path: str) -> Dict[str, Set[str]]:
    """Map cluster id -> set of genome names based on Ndb primary_cluster column."""
    clusters: Dict[str, Set[str]] = defaultdict(set)
    with open(ndb_path, newline="") as handle:
        reader = csv.DictReader(handle)
        if "primary_cluster" not in reader.fieldnames:
            raise ValueError("Ndb file must contain a 'primary_cluster' column.")
        for row in reader:
            genome = row.get("reference") or row.get("querry") or row.get("query")
            cluster = row["primary_cluster"]
            if genome:
                clusters[cluster].add(os.path.basename(genome))
    return clusters


def load_genome_metrics(genome_info_path: str) -> Dict[str, Dict[str, float]]:
    """Return genome -> metrics (completeness, contamination) from genomeInfo.csv."""
    metrics: Dict[str, Dict[str, float]] = {}
    with open(genome_info_path, newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"genome", "completeness", "contamination"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("genomeInfo file must contain 'genome', 'completeness', and 'contamination' columns.")
        for row in reader:
            genome = os.path.basename(row["genome"])
            try:
                comp = float(row["completeness"])
                cont = float(row["contamination"])
            except (TypeError, ValueError):
                continue
            metrics[genome] = {"completeness": comp, "contamination": cont}
    return metrics


def compute_binning_cluster_counts(clusters: Dict[str, Set[str]]) -> Dict[str, int]:
    """Return counts of clusters containing at least one genome per binning type."""
    counts: Dict[str, Set[str]] = defaultdict(set)
    for cluster_id, genomes in clusters.items():
        for genome in genomes:
            counts[binning_type(genome)].add(cluster_id)
    return {bt: len(ids) for bt, ids in counts.items()}


def _binning_sort_key(bt: str):
    """Order: SIN first, then R2..Rn, then ALL, then others alphabetically."""
    parts = bt.split("_")
    tag = parts[1] if len(parts) > 1 else parts[0]
    if tag == "SIN":
        return (0, 0, bt)
    match = re.fullmatch(r"R(\d+)", tag)
    if match:
        return (1, int(match.group(1)), bt)
    if tag == "ALL":
        return (2, 0, bt)
    return (3, 0, bt)


def compute_binning_metrics(
    clusters: Dict[str, Set[str]], metrics: Dict[str, Dict[str, float]]
) -> Dict[str, Dict[str, float]]:
    """
    Compute per-binning-type statistics:
    - clusters: number of clusters containing that binning type
    - completeness: mean of per-cluster mean completeness values
    - contamination: mean of per-cluster mean contamination values
    """
    completeness_cluster_means: Dict[str, list] = defaultdict(list)
    contamination_cluster_means: Dict[str, list] = defaultdict(list)

    for cluster_id, genomes in clusters.items():
        per_type_comp: Dict[str, list] = defaultdict(list)
        per_type_cont: Dict[str, list] = defaultdict(list)
        for genome in genomes:
            if genome not in metrics:
                continue
            bt = binning_type(genome)
            per_type_comp[bt].append(metrics[genome]["completeness"])
            per_type_cont[bt].append(metrics[genome]["contamination"])
        for bt, vals in per_type_comp.items():
            completeness_cluster_means[bt].append(mean(vals))
        for bt, vals in per_type_cont.items():
            contamination_cluster_means[bt].append(mean(vals))

    stats: Dict[str, Dict[str, float]] = {}
    for bt in set(list(completeness_cluster_means.keys()) + list(contamination_cluster_means.keys())):
        cvals = completeness_cluster_means.get(bt, [])
        contvals = contamination_cluster_means.get(bt, [])
        stats[bt] = {
            "clusters": len(cvals),
            "completeness": mean(cvals) if cvals else float("nan"),
            "contamination": mean(contvals) if contvals else float("nan"),
        }
    return stats


def write_binning_metrics(stats: Dict[str, Dict[str, float]], output: str) -> None:
    fieldnames = ["binning_type", "clusters", "completeness", "contamination"]
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for bt in sorted(stats, key=_binning_sort_key):
            row = {"binning_type": bt}
            row.update(stats[bt])
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    metrics_out = os.path.join(args.output_dir, "binning_metrics.csv")
    clusters = load_cluster_members(args.ndb)
    metrics = load_genome_metrics(args.genome_info)
    stats = compute_binning_metrics(clusters, metrics)
    write_binning_metrics(stats, metrics_out)
    print(f"Wrote binning metrics for {len(stats)} binning types to {metrics_out}")


if __name__ == "__main__":
    main()
