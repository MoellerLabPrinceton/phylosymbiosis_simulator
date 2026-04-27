#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def build_output_prefix(
    root_richness: int,
    max_community_size: int,
    n_species_pool: int,
    geom_p: float,
    brownian_sigma: float,
    extinction_rate: float,
    colonization_rate: float,
    n_reps: int,
    n_permutations: int,
    hc_method: str,
    random_seed,
) -> str:
    seed_str = "none" if random_seed is None else str(random_seed)
    return (
        f"root{root_richness}"
        f"_max{max_community_size}"
        f"_pool{n_species_pool}"
        f"_geomp{geom_p}"
        f"_sigma{brownian_sigma}"
        f"_ext{extinction_rate}"
        f"_col{colonization_rate}"
        f"_reps{n_reps}"
        f"_perm{n_permutations}"
        f"_hc{hc_method}"
        f"_seed{seed_str}"
    )


def make_summary_line_plots(summary_files, output_dir: Path) -> None:
    records = []

    pattern = re.compile(
        r"root(?P<root>\d+)_max(?P<max>\d+)_pool(?P<pool>\d+)_geomp(?P<geomp>[\d.]+)"
        r"_sigma(?P<sigma>[\d.]+)_ext(?P<ext>[\d.]+)_col(?P<col>[\d.]+)_reps(?P<reps>\d+)"
        r"_perm(?P<perm>\d+)_hc(?P<hc>[A-Za-z]+)_seed(?P<seed>.+?)\.tsv$"
    )

    for path in summary_files:
        match = pattern.search(path.name)
        if not match:
            print(f"Skipping summary file with unexpected name: {path}")
            continue

        meta = match.groupdict()
        df = pd.read_csv(path, sep="\t")
        df["root_richness"] = int(meta["root"])
        df["max_community_size"] = int(meta["max"])
        df["sigma"] = float(meta["sigma"])
        df["source_file"] = path.name
        records.append(df)

    if not records:
        print("No valid summary files found for summary line plotting.")
        return

    summary = pd.concat(records, ignore_index=True)

    method_order = [
        "Mantel on raw Bray-Curtis",
        "Mantel on cophenetic Bray-Curtis clustering",
        "Normalized Robinson-Foulds",
    ]

    label_map = {
        "Mantel on raw Bray-Curtis": "Mantel raw",
        "Mantel on cophenetic Bray-Curtis clustering": "Mantel cophenetic",
        "Normalized Robinson-Foulds": "Normalized RF",
    }

    color_map = {
        "Mantel on raw Bray-Curtis": "tab:blue",
        "Mantel on cophenetic Bray-Curtis clustering": "tab:orange",
        "Normalized Robinson-Foulds": "tab:green",
    }

    plot_data_path = output_dir / "summary_line_plot_data.tsv"
    summary[
        [
            "sigma",
            "root_richness",
            "max_community_size",
            "test",
            "positive_count",
            "null_count",
            "negative_count",
            "mean_statistic",
            "mean_p_value",
            "source_file",
        ]
    ].sort_values(["sigma", "max_community_size", "test"]).to_csv(
        plot_data_path, sep="\t", index=False
    )

    for sigma in sorted(summary["sigma"].unique()):
        sub = summary[summary["sigma"] == sigma].copy()
        sub = sub.sort_values(["max_community_size", "test"])

        fig, ax = plt.subplots(figsize=(7, 5))
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        for method in method_order:
            m = sub[sub["test"] == method].sort_values("max_community_size")
            if m.empty:
                continue

            x = m["max_community_size"].to_numpy()
            y_pos = m["positive_count"].to_numpy()
            y_null = m["null_count"].to_numpy()

            ax.plot(
                x,
                y_pos,
                linestyle="-",
                linewidth=1.8,
                color=color_map[method],
                label=f"{label_map[method]} positive",
            )
            ax.plot(
                x,
                y_null,
                linestyle="--",
                linewidth=1.8,
                color=color_map[method],
                label=f"{label_map[method]} null",
            )

        ax.set_xlabel("Maximum community size")
        ax.set_ylabel("Count across simulations")
        ax.set_ylim(0, max(100, int(sub[["positive_count", "null_count"]].to_numpy().max())))

        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(1.0)

        ax.tick_params(axis="both", colors="black", width=1.0)
        ax.grid(False)

        legend = ax.legend(
            loc="lower right",
            frameon=True,
            fancybox=False,
            edgecolor="black",
            fontsize=8,
        )
        legend.get_frame().set_linewidth(1.0)

        fig.tight_layout()

        sigma_str = str(sigma).replace(".", "p")
        out_path = output_dir / f"summary_line_plot_sigma{sigma_str}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches="tight")
        plt.close(fig)

        print(f"Saved summary line plot to: {out_path}")

    print(f"Saved summary line plot data to: {plot_data_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper to run the sigfig simulation script for multiple root-richness values. "
            "Each run writes results, summary, and a percent-significant PDF figure. "
            "After all runs, the wrapper also writes summary line plots."
        )
    )
    parser.add_argument(
        "--sim-script",
        default="simulate_microbiome_phylosymbiosis_main.py",
        help="Path to the simulation script",
    )
    parser.add_argument(
        "--tree",
        required=True,
        help="Path to Newick tree",
    )
    parser.add_argument(
        "--root-richnesses",
        nargs="+",
        type=int,
        default=[10, 50, 100, 500],
        help="Root richness values to test",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of CPU cores to pass to each simulation run",
    )
    parser.add_argument(
        "--n-reps",
        type=int,
        default=100,
        help="Number of simulation replicates per parameter set",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=999,
        help="Number of permutations per replicate",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=1,
        help="Random seed passed to each simulation run",
    )
    parser.add_argument(
        "--n-species-pool",
        type=int,
        default=1000,
        help="Total species pool size",
    )
    parser.add_argument(
        "--geom-p",
        type=float,
        default=0.2,
        help="Geometric distribution p for initial abundances",
    )
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
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for output files",
    )
    parser.add_argument(
        "--skip-summary-line-plots",
        action="store_true",
        help="Do not make the aggregate summary line plots after the simulation runs finish",
    )
    args = parser.parse_args()

    sim_script = Path(args.sim_script)
    tree_path = Path(args.tree)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not sim_script.exists():
        raise FileNotFoundError(f"Simulation script not found: {sim_script}")
    if not tree_path.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_path}")

    summary_files = []

    for root_richness in args.root_richnesses:
        max_community_size = int(round(root_richness * 1.5))

        prefix = build_output_prefix(
            root_richness=root_richness,
            max_community_size=max_community_size,
            n_species_pool=args.n_species_pool,
            geom_p=args.geom_p,
            brownian_sigma=args.brownian_sigma,
            extinction_rate=args.extinction_rate,
            colonization_rate=args.colonization_rate,
            n_reps=args.n_reps,
            n_permutations=args.n_permutations,
            hc_method=args.hc_method,
            random_seed=args.random_seed,
        )

        results_out = output_dir / f"simulation_results_{prefix}.tsv"
        summary_out = output_dir / f"simulation_summary_{prefix}.tsv"
        figure_out = output_dir / f"simulation_summary_{prefix}.percent_significant.pdf"

        cmd = [
            sys.executable,
            str(sim_script),
            "--tree", str(tree_path),
            "--n-reps", str(args.n_reps),
            "--n-jobs", str(args.n_jobs),
            "--n-permutations", str(args.n_permutations),
            "--random-seed", str(args.random_seed),
            "--n-species-pool", str(args.n_species_pool),
            "--root-richness", str(root_richness),
            "--max-community-size", str(max_community_size),
            "--geom-p", str(args.geom_p),
            "--brownian-sigma", str(args.brownian_sigma),
            "--extinction-rate", str(args.extinction_rate),
            "--colonization-rate", str(args.colonization_rate),
            "--hc-method", str(args.hc_method),
            "--results-out", str(results_out),
            "--summary-out", str(summary_out),
            "--figure-out", str(figure_out),
        ]

        print("\nRunning:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)

        summary_files.append(summary_out)

    if not args.skip_summary_line_plots:
        make_summary_line_plots(summary_files, output_dir)

    print("\nAll runs completed.")


if __name__ == "__main__":
    main()
