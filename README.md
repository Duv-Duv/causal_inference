# Network causal inference — consolidated R simulations

This directory contains the **two source files** needed for the results simulations and the variance illustration. It can be placed anywhere; it does not require an existing `code/`, `results/`, or `figures/` structure.

## Files

- `network_causal_results.R`: experiments 1–5, academic plots, and optional variance illustration.
- `README.md`: instructions and experiment map.

## Installation

R (recommended >= 4.1) and packages:

```r
install.packages(c("ggplot2", "igraph", "Matrix", "RSpectra"))
```

`Matrix` is usually bundled with R. `igraph`, `Matrix`, and `RSpectra` are required for Experiment 4. The variance illustration alone needs only base R.

## Commands

Run from **any working directory**. By default, outputs are written relative to the current working directory:

```bash
Rscript network_causal_results.R --quick
Rscript network_causal_results.R
Rscript network_causal_results.R --variance
Rscript network_causal_results.R --only-variance
```

To explicitly choose an output directory (created automatically if absent):

```bash
Rscript network_causal_results.R --quick ./my_output
```

`--quick` uses 30 Monte Carlo repetitions for Experiments 1, 2, 3, 5 and 12 for Experiment 4. The full run uses 180 for all experiments. For the variance illustration, `--quick` uses 300 randomizations rather than 6000. `--variance` runs all five experiments **and** the illustration; `--only-variance` runs only the illustration.

**Note:** Experiment 4 with n=10000 and 180 repetitions can be computationally demanding, even with sparse matrices. The quick run is a smoke test, not publication-quality Monte Carlo evidence.

## Generated folders

The script creates these folders when needed:

```text
<output-root>/
  results/
    positivity_draws.csv
    positivity_summary.csv
    positivity_exact_probabilities.csv
    positivity_diagnostic.csv
    density_estimators_summary.csv
    pc_balancing_summary.csv
    scaling_summary.csv
    network_error_summary.csv
  figures/
    exact_exposure_probabilities.pdf/.png
    exact_positivity_rmse.pdf/.png
    direct_density_rmse.pdf/.png
    indirect_unbiased_rmse.pdf/.png
    pc_rank_rmse.pdf/.png
    sparsity_scaling.pdf/.png
    network_measurement_error.pdf/.png
    academic/  # only when --variance or --only-variance
      indirect_variance_explosion.pdf/.png
```

## Experiments and preserved revisions

| Experiment | Focus | Version retained |
|---|---|---|
| 1 | Exact exposure probabilities, positivity, HT vs unweighted | Latest seed-25 revision: n=1000, K=100, ten densities, common-support diagnostics |
| 2 | Direct HT/Hajek and natural indirect estimator | Latest main suite |
| 3 | PC balancing, ranks 1/2/4 | Latest main suite; eigen decomposition reused across ranks |
| 4 | Sparse-sequence scaling | Latest original-style version, n = 100, 1000, 10000; sparse igraph/Matrix/RSpectra |
| 5 | Network measurement error | Latest main suite |
| Illustration | Variance explosion of natural indirect estimator | Original standalone base-R figure integrated as an optional mode |

