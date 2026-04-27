# Microbiome Phylosymbiosis Simulation

This repository contains two scripts:

1. A simulation engine that evolves microbiome communities along a host phylogeny and tests for phylosymbiosis  
2. A wrapper that runs the simulation across multiple root richness values and organizes outputs  

---

## Scripts

### 1. simulate_microbiome_phylosymbiosis_with_sigfig.py

Main simulation script.

Simulates microbiome evolution along a host tree, computes dissimilarities, and evaluates phylosymbiosis using three approaches:

- Mantel test on raw Bray Curtis distances  
- Mantel test on cophenetic distances from hierarchical clustering  
- Normalized Robinson Foulds comparison between microbial and host trees  

Outputs:

- Per replicate results table  
- Summary table across replicates  
- PDF figure showing percent significant results  

---

### 2. run_root_richness_wrapper_sigfig.py

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
