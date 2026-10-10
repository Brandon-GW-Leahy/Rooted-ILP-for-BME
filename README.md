# Balanced Minimum Evolution: ILP Comparison

This repository contains `ilp_comparison.py`, a Python implementation for comparing two integer linear programming (ILP) formulations of the **Balanced Minimum Evolution (BME)** problem for phylogenetic tree reconstruction:

1. **Rooted-tree ILP:** The compact binary-encoding formulation developed by Brandon Leahy and Pavel Skums, in which taxa are represented by binary root-to-leaf paths and attachment depths.
2. **Path-length ILP:** The formulation based on the pairwise path-length matrix of an unrooted binary tree, following Catanzaro et al., *Optimizing over path-length matrices of unrooted binary trees* (2025).

Both formulations are solved with Gurobi on the same dissimilarity matrix in each trial. The script records the best objective value found and the solver runtime for each formulation.

For details of the rooted formulation and the experimental comparisons, see Leahy and Skums, *A Compact Rooted-Tree ILP Formulation for Balanced Minimum Evolution Problem* (October 2026).

## Requirements

- **Python 3.11** (used for the experiments reported in the paper).
- **Gurobi Optimizer** and an appropriate Gurobi license. The paper's experiments used Gurobi 10.0.2. A license that supports the intended model sizes is necessary for larger instances.
- Python packages: `gurobipy`, `numpy`, and `openpyxl`.

Install the Python dependencies with:

```bash
python -m pip install gurobipy numpy openpyxl
```

Gurobi must also be correctly licensed in the environment where the script is run. The `numpy` package is imported by the script, although it is not used by the current experiment loop.

## Running the script

From the directory containing `ilp_comparison.py`, run:

```bash
python ilp_comparison.py
```

The script will ask for four values:

| Prompt | Description |
| --- | --- |
| `Number of taxa` | Number of taxa, `n`, in each input matrix. |
| `Maximum length of rooted tree` | Maximum rooted-tree depth, `L`, used by the rooted ILP. |
| `Number of trials` | Number of independent input matrices to generate and compare. |
| `Time limit (hrs)` | Integer number of hours allowed **per formulation, per trial**. |

For example, a small test can use `n = 6`, `L = 3`, `1` trial, and a time limit of `1` hour. The path-length formulation automatically uses the bound `R = 2L - 1`; it is not requested separately.

**Note:** The two models run sequentially, each with its own time limit. A run of `k` trials can therefore require up to approximately `2 × k × (time limit)` hours, plus model-building and file-writing overhead. Larger instances may require considerable memory and computation time. Choose `L` large enough to allow a feasible rooted tree; in particular, `n` cannot exceed `2^L`.

## Selecting the input data

The script expects a symmetric, `n × n` dissimilarity matrix `D` with a zero diagonal. Both models receive the **same** `D` in each trial.

### 1. Random Euclidean-distance matrices (default)

As uploaded, the main trial loop uses:

```python
D = build_distance_matrix(n)
```

This samples `n` points uniformly from the three-dimensional unit ball and computes their pairwise Euclidean distances. Each trial generates a new matrix.

### 2. Noisy-additive matrices

The function `generate_random_dissimilarity_matrix(n)` is also provided. It generates a random binary phylogenetic tree, calculates the associated unrooted leaf-to-leaf path distances, and adds independent integer noise sampled from `0` through `9` to each off-diagonal pair.

To use this generator, replace the matrix-generation line **inside the main trial loop** with:

```python
D = generate_random_dissimilarity_matrix(n)
```

### 3. Custom or published benchmark matrices

To run a particular matrix, replace the matrix-generation line **inside the main trial loop** with a manually specified Python list of lists. For example, the following is an *illustrative* four-taxon matrix, **not** one of the published benchmarks:

```python
D = [
    [0, 2, 3, 3],
    [2, 0, 3, 3],
    [3, 3, 0, 2],
    [3, 3, 2, 0],
]
```

Set `Number of taxa` to `4` for this example, and use `Number of trials = 1` when testing a single fixed matrix. To evaluate several distinct supplied matrices, modify the loop to select a different `D` on each iteration. Ensure the entered `n` matches the matrix dimensions.

**Benchmark-data availability:** The RIM (random integer matrix), RDSM (random doubly stochastic matrix), and Phylobench instances discussed in the paper were manually entered into the script for the corresponding experiments. **Those benchmark matrices are not included in this repository.** The uploaded version does not contain a benchmark-data directory, a benchmark-selection option, or a file-loading interface. Reproducing the paper's benchmark runs therefore requires obtaining the appropriate original matrices separately and inserting them into the script (or adding your own data loader). Running the unmodified script will generate new Euclidean instances; it will **not** reproduce the RIM, RDSM, or Phylobench tables.

The RIM and RDSM benchmark families are discussed in Catanzaro, Dehaybe, and Pesenti, *New heuristics for phylogeny estimation under the balanced minimum evolution criterion* (2025). See the accompanying Leahy–Skums paper for the specific benchmark instance names and Phylobench experiments.

## Results

After all trials complete, the script writes an Excel workbook in the current working directory:

```text
ILP_comparison_results.xlsx
```

The `Results` worksheet contains one row per trial and four columns:

| Column | Meaning |
| --- | --- |
| `Best Rooted ILP Objective` | Best BME objective value found by the rooted-tree model. |
| `Rooted ILP Run Time` | Gurobi runtime for the rooted-tree model, in seconds. |
| `Best Path-Length ILP Objective` | Best BME objective value found by the path-length model. |
| `Path-Length ILP Run Time` | Gurobi runtime for the path-length model, in seconds. |

If Gurobi finds no feasible solution within the allotted time, the corresponding objective is stored as `None` in Python and normally appears as a blank Excel cell. A reported objective is the **best feasible value found**, not necessarily a proven global optimum; the script does not save solver optimality statuses or optimality gaps.

**The output workbook is overwritten on each run.** Rename or move it if you want to retain results from multiple experiments. The script does not currently export the inferred trees, input matrices, or Gurobi models.

## Experimental settings used in the paper

The synthetic experiments reported in the accompanying paper use the following settings:

| Number of taxa (`n`) | Rooted depth (`L`) | Path-length bound (`R = 2L - 1`) | Time limit per model |
| ---: | ---: | ---: | ---: |
| 20 | 5 | 9 | 1 hour |
| 30 | 6 | 11 | 5 hours |
| 40 | 6 | 11 | 12 hours |

The paper additionally considers Phylobench instances with `n = 30` (5-hour limit) and `n = 45` (14-hour limit), both with `L = 6` and `R = 11`. The published experiments ran on Ubuntu 22.04 LTS using Python 3.11 and Gurobi 10.0.2, on a workstation with two Intel Xeon Platinum 8268 processors and 251 GiB of RAM. Runtimes may differ substantially on other systems.

For reproducible **newly generated** instances, a fixed Python random seed can be added near the beginning of the script, for example:

```python
random.seed(42)
```

This controls the random input generation but does not recover omitted benchmark matrices or guarantee identical optimization runtimes.

## References

- Brandon Leahy and Pavel Skums. *A Compact Rooted-Tree ILP Formulation for Balanced Minimum Evolution Problem*. Manuscript, October 2026.
- D. Catanzaro, R. Pesenti, A. Sapucaia Barboza, and L. Wolsey. *Optimizing over path-length matrices of unrooted binary trees*. Mathematical Programming, 2025.
- D. Catanzaro, H. Dehaybe, and R. Pesenti. *New heuristics for phylogeny estimation under the balanced minimum evolution criterion*. Bioinformatics, 2025.
