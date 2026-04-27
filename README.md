# Microbiome Phylosymbiosis Simulation

This repository contains two scripts:

1. An engine that simulates microbiome composition changes along a host phylogeny and tests for phylosymbiosis  
2. A wrapper that runs the simulation across multiple root richness values and organizes outputs  

---

## Scripts

### 1. simulate_microbiome_phylosymbiosis_main.py

Main simulation script.

Simulates microbiome composition changes along a host tree, computes dissimilarities, and evaluates phylosymbiosis using three approaches:

- Mantel test on raw Bray Curtis distances  
- Mantel test on cophenetic distances from hierarchical clustering  
- Normalized Robinson Foulds comparison between microbial and host trees  

Outputs:

- Per replicate results table  
- Summary table across replicates  
- PDF figure showing percent significant results  

---

### 2. run_root_richness_wrapper_main.py

Wrapper script for batch experiments.

Runs the simulation for multiple root richness values and automatically sets:

- max-community-size = 1.5 × root richness  

Generates parameter encoded filenames and corresponding figures.

---

## Installation

Requires Python 3.8 or higher.

Install dependencies:

```bash
pip install numpy scipy biopython matplotlib
```

---

## Input

A Newick formatted phylogenetic tree with named tips.

Example:

```bash
--tree new_tree.nwk
```

---

## Usage

### Running the simulation script

Basic example:

```bash
python simulate_microbiome_phylosymbiosis_main.py \
  --tree new_tree.nwk \
  --n-reps 100 \
  --n-jobs 8 \
  --root-richness 100 \
  --max-community-size 150 \
  --results-out results.tsv \
  --summary-out summary.tsv
```

Outputs:

- results.tsv  
- summary.tsv  
- summary.percent_significant.pdf  

Optional:

```bash
--figure-out custom_figure.pdf
```

---

### Running the wrapper script

```bash
python run_root_richness_wrapper_main.py \
  --tree new_tree.nwk \
  --n-reps 100 \
  --n-jobs 8
```

Default root richness values:

- 10  
- 50  
- 100  
- 500  

Each run automatically sets:

- max-community-size = 1.5 × root richness  

Example output filenames:

```
simulation_results_root100_max150_pool1000_geomp0.2_sigma0.35_ext0.08_col0.04_reps100_perm999_hcaverage_seed1.tsv

simulation_summary_root100_max150_pool1000_geomp0.2_sigma0.35_ext0.08_col0.04_reps100_perm999_hcaverage_seed1.tsv

simulation_summary_root100_max150_pool1000_geomp0.2_sigma0.35_ext0.08_col0.04_reps100_perm999_hcaverage_seed1.percent_significant.pdf
```

---

## Parameters

### Core simulation parameters

- `--tree`  
  Path to Newick tree  

- `--n-reps`  
  Number of replicate simulations  

- `--n-jobs`  
  Number of CPU cores  

- `--n-permutations`  
  Number of permutations for significance tests  

- `--random-seed`  
  Controls reproducibility  

---

### Community initialization

- `--n-species-pool`  
  Total number of possible species  

- `--root-richness`  
  Number of species present at the root  

- `--geom-p`  
  Parameter of geometric distribution for initial abundances  
  Lower values produce more skewed communities  

---

### Evolution parameters

- `--brownian-sigma`  
  Strength of stochastic change in log abundance  

- `--extinction-rate`  
  Rate at which taxa are lost  

- `--colonization-rate`  
  Rate at which new taxa appear  

---

### Community constraint

- `--max-community-size`  
  Maximum number of species allowed per community  

  If exceeded:
  - taxa are probabilistically removed  
  - removal probability is proportional to abundance  

---

### Clustering

- `--hc-method`  

Options include:

- average  
- complete  
- single  
- ward  

---

### Output control

- `--results-out`  
  Per replicate results file  

- `--summary-out`  
  Summary file across replicates  

- `--figure-out`  
  Output PDF for percent significant plot  

---

## Model overview

### Initialization

- Root community is sampled from a species pool  
- Abundances drawn from a geometric distribution  

---

### Evolution along branches

For each branch:

- Abundances evolve via Brownian motion on log scale  
- Taxa go extinct with probability determined by extinction rate  
- New taxa colonize with probability determined by colonization rate  
- Community size is constrained if needed  

---

## Phylosymbiosis tests

### Mantel tests

- Compare host phylogenetic distances to microbiome distances  
- Applied to:
  - raw Bray Curtis  
  - cophenetic distances  

---

### Normalized Robinson Foulds

- Compares microbial clustering topology to host tree  
- Uses permutation testing  
- Lower values indicate stronger similarity  

---

## Output interpretation

Summary table reports:

- Mean and median test statistics  
- Mean p values  
- Counts of:
  - positive results  
  - null results  
  - negative results  

---

## Percent significant figure

The PDF shows:

- Fraction of simulations with significant results  
- Separate curves for positive and negative signal  

---

## Notes

- Reproducibility depends on the random seed  
- Higher brownian-sigma increases divergence among communities  
- Lower geom-p produces stronger dominance patterns  
- Larger max-community-size reduces pruning effects  
