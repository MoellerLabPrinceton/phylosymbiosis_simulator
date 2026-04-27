#!/usr/bin/env python3

import argparse
import csv
import math
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List, Tuple, Optional
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class SimulationParams:
    tree_path: str
    n_species_pool: int = 1000
    root_richness: int = 100
    max_community_size: Optional[int] = None
    geom_p: float = 0.2
    brownian_sigma: float = 0.35
    extinction_rate: float = 0.08
    colonization_rate: float = 0.04
    min_positive_abundance: float = 1e-12
    n_reps: int = 100
    n_permutations: int = 999
    random_seed: Optional[int] = 1
    hc_method: str = "average"
    n_jobs: int = 1
    results_out: str = "simulation_results.tsv"
    summary_out: str = "simulation_summary.tsv"
    figure_out: Optional[str] = None


def normalize_profile(profile: np.ndarray) -> np.ndarray:
    total = profile.sum()
    if total <= 0:
        raise ValueError("Profile has zero total abundance.")
    return profile / total


def bray_curtis(u: np.ndarray, v: np.ndarray) -> float:
    denom = u.sum() + v.sum()
    if denom == 0:
        return 0.0
    return np.abs(u - v).sum() / denom


def pairwise_bray_curtis(profiles: np.ndarray) -> np.ndarray:
    n = profiles.shape[0]
    out = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = bray_curtis(profiles[i], profiles[j])
            out[i, j] = d
            out[j, i] = d
    return out


def lower_triangle_values(mat: np.ndarray) -> np.ndarray:
    idx = np.tril_indices_from(mat, k=-1)
    return mat[idx]


def mantel_test(
    mat1: np.ndarray,
    mat2: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
) -> Tuple[float, float]:
    if mat1.shape != mat2.shape:
        raise ValueError("Matrices must have same shape.")
    if mat1.shape[0] != mat1.shape[1]:
        raise ValueError("Mantel test requires square matrices.")

    x = lower_triangle_values(mat1)
    y = lower_triangle_values(mat2)

    if len(x) == 0 or len(y) == 0:
        return 0.0, 1.0

    if x.std(ddof=1) == 0 or y.std(ddof=1) == 0:
        return 0.0, 1.0

    observed = np.corrcoef(x, y)[0, 1]
    n = mat1.shape[0]
    perm_stats = np.empty(n_permutations, dtype=float)

    for k in range(n_permutations):
        perm = rng.permutation(n)
        y_perm = lower_triangle_values(mat2[perm][:, perm])

        if len(y_perm) == 0 or y_perm.std(ddof=1) == 0:
            perm_stats[k] = 0.0
        else:
            perm_stats[k] = np.corrcoef(x, y_perm)[0, 1]

    p = (np.sum(np.abs(perm_stats) >= abs(observed)) + 1) / (len(perm_stats) + 1)
    return observed, p


def classify_mantel_evidence(r: float, p: float, alpha: float = 0.05) -> str:
    if np.isnan(r) or np.isnan(p):
        return "null"
    if p < alpha and r > 0:
        return "positive"
    if p < alpha and r < 0:
        return "negative"
    return "null"


def classify_rf_evidence(
    obs_nrf: float,
    p_lower: float,
    p_upper: float,
    alpha: float = 0.05,
) -> str:
    if np.isnan(obs_nrf) or np.isnan(p_lower) or np.isnan(p_upper):
        return "null"
    if p_lower < alpha:
        return "positive"
    if p_upper < alpha:
        return "negative"
    return "null"


def get_patristic_distance_matrix(tree) -> Tuple[List[str], np.ndarray]:
    tips = tree.get_terminals()
    names = [tip.name for tip in tips]
    n = len(tips)
    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            d = tree.distance(tips[i], tips[j])
            dist[i, j] = d
            dist[j, i] = d

    return names, dist


def enforce_max_richness_weighted(
    profile: np.ndarray,
    max_k: Optional[int],
    rng: np.random.Generator,
) -> np.ndarray:
    if max_k is None:
        return profile

    present_idx = np.where(profile > 0)[0]
    richness = len(present_idx)

    if richness <= max_k:
        return profile

    abund = profile[present_idx].astype(float)
    total = abund.sum()

    if total <= 0:
        keep_idx = rng.choice(present_idx, size=max_k, replace=False)
    else:
        probs = abund / total
        keep_idx = rng.choice(present_idx, size=max_k, replace=False, p=probs)

    new_profile = np.zeros_like(profile)
    new_profile[keep_idx] = profile[keep_idx]
    return new_profile


def initialize_root_profile(params: SimulationParams, rng: np.random.Generator) -> np.ndarray:
    profile = np.zeros(params.n_species_pool, dtype=float)
    chosen = rng.choice(params.n_species_pool, size=params.root_richness, replace=False)
    abundances = rng.geometric(p=params.geom_p, size=params.root_richness).astype(float)
    profile[chosen] = abundances

    if params.max_community_size is not None:
        profile = enforce_max_richness_weighted(profile, params.max_community_size, rng)

    return normalize_profile(profile)


def evolve_profile_along_branch(
    parent_profile: np.ndarray,
    branch_length: float,
    params: SimulationParams,
    rng: np.random.Generator,
) -> np.ndarray:
    if branch_length is None:
        branch_length = 1.0
    if branch_length < 0:
        raise ValueError("Negative branch length detected.")

    profile = parent_profile.copy()
    present = profile > 0
    absent = ~present

    if present.any():
        sigma = params.brownian_sigma * math.sqrt(max(branch_length, 0.0))
        log_abund = np.log(profile[present] + params.min_positive_abundance)
        log_abund += rng.normal(0.0, sigma, size=log_abund.shape[0])
        updated = np.exp(log_abund)

        p_ext = 1.0 - math.exp(-params.extinction_rate * branch_length)
        extinct = rng.random(updated.shape[0]) < p_ext
        updated[extinct] = 0.0
        profile[present] = updated

    absent_indices = np.where(absent)[0]
    if len(absent_indices) > 0:
        p_col = 1.0 - math.exp(-params.colonization_rate * branch_length)
        colonize = rng.random(len(absent_indices)) < p_col
        chosen_absent = absent_indices[colonize]
        if len(chosen_absent) > 0:
            new_abund = rng.geometric(
                p=params.geom_p, size=len(chosen_absent)
            ).astype(float) / 10.0
            profile[chosen_absent] = new_abund

    profile[profile < params.min_positive_abundance] = 0.0

    if params.max_community_size is not None:
        profile = enforce_max_richness_weighted(
            profile,
            params.max_community_size,
            rng,
        )

    if profile.sum() <= 0:
        fallback = rng.choice(params.n_species_pool, size=1, replace=False)
        profile[fallback] = 1.0

    return normalize_profile(profile)


def simulate_tip_profiles(
    tree,
    params: SimulationParams,
    rng: np.random.Generator,
) -> Dict[str, np.ndarray]:
    root_profile = initialize_root_profile(params, rng)
    tip_profiles: Dict[str, np.ndarray] = {}

    def recurse(clade, parent_profile: np.ndarray) -> None:
        branch_length = clade.branch_length if clade.branch_length is not None else 0.0
        current_profile = evolve_profile_along_branch(parent_profile, branch_length, params, rng)

        if clade.is_terminal():
            if clade.name is None:
                raise ValueError("All terminal tips must have names.")
            tip_profiles[clade.name] = current_profile
            return

        for child in clade.clades:
            recurse(child, current_profile)

    for child in tree.root.clades:
        recurse(child, root_profile)

    return tip_profiles


def cophenetic_from_bray_curtis(bray_mat: np.ndarray, method: str) -> np.ndarray:
    condensed = squareform(bray_mat, checks=False)
    Z = linkage(condensed, method=method)
    n = bray_mat.shape[0]
    coph = np.zeros((n, n), dtype=float)

    clusters = {i: {i} for i in range(n)}

    for row_idx, row in enumerate(Z):
        a = int(row[0])
        b = int(row[1])
        dist = float(row[2])
        new_id = n + row_idx

        left = clusters[a]
        right = clusters[b]
        members = left | right

        for i in left:
            for j in right:
                coph[i, j] = dist
                coph[j, i] = dist

        clusters[new_id] = members

    return coph


def canonical_split(subset: set, all_leaves: set) -> frozenset:
    other = all_leaves - subset
    subset_sorted = tuple(sorted(subset))
    other_sorted = tuple(sorted(other))

    if len(subset) < len(other):
        return frozenset(subset)
    if len(other) < len(subset):
        return frozenset(other)
    return frozenset(subset if subset_sorted <= other_sorted else other)


def get_host_splits(tree, tip_names: List[str]) -> set:
    all_leaves = set(tip_names)
    splits = set()

    def recurse(clade):
        if clade.is_terminal():
            return {clade.name}
        desc = set()
        for child in clade.clades:
            desc |= recurse(child)
        if 1 < len(desc) < len(all_leaves) - 1:
            splits.add(canonical_split(desc, all_leaves))
        return desc

    recurse(tree.root)
    return splits


def linkage_splits(Z: np.ndarray, tip_names: List[str]) -> set:
    n = len(tip_names)
    all_leaves = set(tip_names)
    clusters = {i: {tip_names[i]} for i in range(n)}
    splits = set()

    for row_idx, row in enumerate(Z):
        a = int(row[0])
        b = int(row[1])
        new_id = n + row_idx
        members = clusters[a] | clusters[b]
        if 1 < len(members) < len(all_leaves) - 1:
            splits.add(canonical_split(members, all_leaves))
        clusters[new_id] = members

    return splits


def normalized_robinson_foulds(splits1: set, splits2: set) -> float:
    denom = len(splits1) + len(splits2)
    if denom == 0:
        return 1.0
    rf = len(splits1 - splits2) + len(splits2 - splits1)
    return rf / denom


def rf_permutation_test(
    microbial_splits: set,
    host_splits: set,
    host_tip_names: List[str],
    n_permutations: int,
    rng: np.random.Generator,
) -> Tuple[float, float, float]:
    obs = normalized_robinson_foulds(microbial_splits, host_splits)

    if np.isnan(obs):
        return 1.0, 1.0, 1.0

    perm_stats = np.empty(n_permutations, dtype=float)
    original_names = list(host_tip_names)
    all_leaves = set(original_names)

    for k in range(n_permutations):
        permuted_names = list(rng.permutation(original_names))
        mapping = {old: new for old, new in zip(original_names, permuted_names)}

        perm_splits = set()
        for split in host_splits:
            renamed = {mapping[x] for x in split}
            perm_splits.add(canonical_split(renamed, all_leaves))

        perm_stats[k] = normalized_robinson_foulds(microbial_splits, perm_splits)

    p_lower = (np.sum(perm_stats <= obs) + 1) / (n_permutations + 1)
    p_upper = (np.sum(perm_stats >= obs) + 1) / (n_permutations + 1)

    return obs, p_lower, p_upper


def write_tsv(path: str, fieldnames: List[str], rows: List[dict]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_results(rows: List[dict]) -> List[dict]:
    out = []

    def mean(vals):
        vals = [x for x in vals if not np.isnan(x)]
        return float(np.mean(vals)) if vals else np.nan

    def median(vals):
        vals = [x for x in vals if not np.isnan(x)]
        return float(np.median(vals)) if vals else np.nan

    raw_r = [r["r_raw"] for r in rows]
    raw_p = [r["p_raw"] for r in rows]
    raw_e = [r["evidence_raw"] for r in rows]

    out.append({
        "test": "Mantel on raw Bray-Curtis",
        "mean_statistic": mean(raw_r),
        "median_statistic": median(raw_r),
        "mean_p_value": mean(raw_p),
        "positive_count": raw_e.count("positive"),
        "null_count": raw_e.count("null"),
        "negative_count": raw_e.count("negative"),
    })

    coph_r = [r["r_cophenetic"] for r in rows]
    coph_p = [r["p_cophenetic"] for r in rows]
    coph_e = [r["evidence_cophenetic"] for r in rows]

    out.append({
        "test": "Mantel on cophenetic Bray-Curtis clustering",
        "mean_statistic": mean(coph_r),
        "median_statistic": median(coph_r),
        "mean_p_value": mean(coph_p),
        "positive_count": coph_e.count("positive"),
        "null_count": coph_e.count("null"),
        "negative_count": coph_e.count("negative"),
    })

    rf_d = [r["normalized_rf"] for r in rows]
    rf_p = [r["p_rf_lower"] for r in rows]
    rf_e = [r["evidence_rf"] for r in rows]

    out.append({
        "test": "Normalized Robinson-Foulds",
        "mean_statistic": mean(rf_d),
        "median_statistic": median(rf_d),
        "mean_p_value": mean(rf_p),
        "positive_count": rf_e.count("positive"),
        "null_count": rf_e.count("null"),
        "negative_count": rf_e.count("negative"),
    })

    return out



def default_figure_path(summary_out: str) -> str:
    path = Path(summary_out)
    if path.suffix:
        return str(path.with_suffix(".percent_significant.pdf"))
    return str(path) + ".percent_significant.pdf"


def plot_percent_significant(summary_rows: List[dict], figure_out: str) -> None:
    methods = [row["test"] for row in summary_rows]
    positive = []
    negative = []

    for row in summary_rows:
        total = row["positive_count"] + row["null_count"] + row["negative_count"]
        if total <= 0:
            positive.append(0.0)
            negative.append(0.0)
        else:
            positive.append(100.0 * row["positive_count"] / total)
            negative.append(100.0 * row["negative_count"] / total)

    labels = [
        "Raw Mantel" if m == "Mantel on raw Bray-Curtis" else
        "Cophenetic Mantel" if m == "Mantel on cophenetic Bray-Curtis clustering" else
        "Normalized RF" if m == "Normalized Robinson-Foulds" else
        m
        for m in methods
    ]

    x = np.arange(len(labels))
    width = 0.36

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.bar(x - width / 2, positive, width, label="Positive")
    ax.bar(x + width / 2, negative, width, label="Negative")

    ax.set_ylabel("Significant simulations (%)")
    ax.set_ylim(0, 100)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")

    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.tick_params(axis="both", colors="black", width=1.0)
    ax.grid(False)

    legend = ax.legend(loc="upper right", frameon=True, fancybox=False, edgecolor="black")
    legend.get_frame().set_linewidth(1.0)

    fig.tight_layout()
    fig.savefig(figure_out, format="pdf", bbox_inches="tight")
    plt.close(fig)

def format_float(x):
    if isinstance(x, float) and np.isnan(x):
        return "nan"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def run_one_simulation(task: Tuple[int, SimulationParams]) -> dict:
    rep_idx, params = task

    if params.random_seed is None:
        rng = np.random.default_rng()
    else:
        rng = np.random.default_rng(params.random_seed + rep_idx)

    tree = Phylo.read(params.tree_path, "newick")
    host_order, host_phylo_dist = get_patristic_distance_matrix(tree)
    host_splits = get_host_splits(tree, host_order)

    tip_profiles = simulate_tip_profiles(tree, params, rng)
    profiles = np.vstack([tip_profiles[name] for name in host_order])

    bray_mat = pairwise_bray_curtis(profiles)
    coph_mat = cophenetic_from_bray_curtis(bray_mat, method=params.hc_method)

    r_raw, p_raw = mantel_test(bray_mat, host_phylo_dist, params.n_permutations, rng)
    r_coph, p_coph = mantel_test(coph_mat, host_phylo_dist, params.n_permutations, rng)

    condensed = squareform(bray_mat, checks=False)
    Z = linkage(condensed, method=params.hc_method)
    microbial_splits = linkage_splits(Z, host_order)

    nrf, p_rf_lower, p_rf_upper = rf_permutation_test(
        microbial_splits,
        host_splits,
        host_order,
        params.n_permutations,
        rng,
    )

    bray_lower = lower_triangle_values(bray_mat)
    phylo_lower = lower_triangle_values(host_phylo_dist)
    coph_lower = lower_triangle_values(coph_mat)

    return {
        "replicate": rep_idx + 1,
        "r_raw": r_raw,
        "p_raw": p_raw,
        "evidence_raw": classify_mantel_evidence(r_raw, p_raw),
        "r_cophenetic": r_coph,
        "p_cophenetic": p_coph,
        "evidence_cophenetic": classify_mantel_evidence(r_coph, p_coph),
        "normalized_rf": nrf,
        "p_rf_lower": p_rf_lower,
        "p_rf_upper": p_rf_upper,
        "evidence_rf": classify_rf_evidence(nrf, p_rf_lower, p_rf_upper),
        "mean_bray_curtis": float(bray_lower.mean()) if len(bray_lower) > 0 else 0.0,
        "mean_richness": float((profiles > 0).sum(axis=1).mean()),
        "host_split_count": len(host_splits),
        "microbial_split_count": len(microbial_splits),
        "bray_variance": float(np.var(bray_lower)) if len(bray_lower) > 0 else 0.0,
        "phylo_variance": float(np.var(phylo_lower)) if len(phylo_lower) > 0 else 0.0,
        "coph_variance": float(np.var(coph_lower)) if len(coph_lower) > 0 else 0.0,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Simulate microbiome evolution on a host tree and test phylosymbiosis."
    )
    parser.add_argument("--tree", required=True, help="Path to Newick tree")
    parser.add_argument("--n-reps", type=int, default=100, help="Number of simulation replicates")
    parser.add_argument("--n-jobs", type=int, default=1, help="Number of CPU cores to use")
    parser.add_argument("--n-permutations", type=int, default=999, help="Permutation count for significance tests")
    parser.add_argument("--random-seed", type=int, default=1, help="Random seed")

    parser.add_argument("--n-species-pool", type=int, default=1000, help="Total species pool size")
    parser.add_argument("--root-richness", type=int, default=100, help="Number of species present at root")
    parser.add_argument(
        "--max-community-size",
        type=int,
        default=None,
        help="Maximum number of species allowed in a community",
    )
    parser.add_argument("--geom-p", type=float, default=0.2, help="Geometric distribution p for initial abundances")
    parser.add_argument(
        "--brownian-sigma",
        type=float,
        default=0.35,
        help="Brownian motion sigma on log abundance per sqrt(branch length)",
    )
    parser.add_argument(
        "--extinction-rate",
        type=float,
        default=0.08,
        help="Extinction rate per unit branch length",
    )
    parser.add_argument(
        "--colonization-rate",
        type=float,
        default=0.04,
        help="Colonization rate per unit branch length",
    )
    parser.add_argument(
        "--hc-method",
        choices=["single", "complete", "average", "weighted", "centroid", "median", "ward"],
        default="average",
        help="Hierarchical clustering method for microbial dendrogram",
    )

    parser.add_argument("--results-out", default="simulation_results.tsv")
    parser.add_argument("--summary-out", default="simulation_summary.tsv")
    parser.add_argument("--figure-out", default=None, help="Path for PDF plot of percent significant results")
    args = parser.parse_args()

    params = SimulationParams(
        tree_path=args.tree,
        n_species_pool=args.n_species_pool,
        root_richness=args.root_richness,
        max_community_size=args.max_community_size,
        geom_p=args.geom_p,
        brownian_sigma=args.brownian_sigma,
        extinction_rate=args.extinction_rate,
        colonization_rate=args.colonization_rate,
        n_reps=args.n_reps,
        n_permutations=args.n_permutations,
        random_seed=args.random_seed,
        hc_method=args.hc_method,
        n_jobs=args.n_jobs,
        results_out=args.results_out,
        summary_out=args.summary_out,
        figure_out=args.figure_out,
    )

    tasks = [(i, params) for i in range(params.n_reps)]

    if params.n_jobs == 1:
        rows = [run_one_simulation(task) for task in tasks]
    else:
        with Pool(processes=params.n_jobs) as pool:
            rows = pool.map(run_one_simulation, tasks)

    rows.sort(key=lambda x: x["replicate"])
    summary_rows = summarize_results(rows)

    result_fields = [
        "replicate",
        "r_raw",
        "p_raw",
        "evidence_raw",
        "r_cophenetic",
        "p_cophenetic",
        "evidence_cophenetic",
        "normalized_rf",
        "p_rf_lower",
        "p_rf_upper",
        "evidence_rf",
        "mean_bray_curtis",
        "mean_richness",
        "host_split_count",
        "microbial_split_count",
        "bray_variance",
        "phylo_variance",
        "coph_variance",
    ]
    summary_fields = [
        "test",
        "mean_statistic",
        "median_statistic",
        "mean_p_value",
        "positive_count",
        "null_count",
        "negative_count",
    ]

    write_tsv(params.results_out, result_fields, rows)
    write_tsv(params.summary_out, summary_fields, summary_rows)

    figure_out = params.figure_out if params.figure_out is not None else default_figure_path(params.summary_out)
    plot_percent_significant(summary_rows, figure_out)

    print("\nPer-replicate results:")
    print("\t".join(result_fields))
    for row in rows:
        print("\t".join(format_float(row[f]) for f in result_fields))

    print("\nSummary across simulations:")
    print("\t".join(summary_fields))
    for row in summary_rows:
        print("\t".join(format_float(row[f]) for f in summary_fields))

    print(f"\nSaved per-replicate results to: {params.results_out}")
    print(f"Saved summary to: {params.summary_out}")
    print(f"Saved percent significant figure to: {figure_out}")


if __name__ == "__main__":
    main()
