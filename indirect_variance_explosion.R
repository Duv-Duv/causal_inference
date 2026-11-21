# Illustrating the variance explosion of the natural indirect-effect estimator
#
# Model:
#   Y_i = b + epsilon_i
#   T_b = b / {n p (1-p)} * sum_j N_j (Z_j - p)
#
# The figure has two panels:
#   (1) The expected conditional variance over an Erdos-Renyi G(n, rho) graph.
#   (2) Assignment distributions of T_b for several graph sizes, conditional
#       on one simulated Erdos-Renyi degree sequence at each size.
#
# Only base R is required.

set.seed(2026)

# Resolve the project root from an optional command-line argument or the
# location of this script, so figures are saved in the project output tree.
script_args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs()
file_arg <- grep("^--file=", all_args, value = TRUE)
project_root <- if (length(script_args) > 0) {
  normalizePath(script_args[1], mustWork = TRUE)
} else if (length(file_arg) > 0) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."),
                mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
figure_dir <- file.path(project_root, "figures", "academic")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Parameters
# -----------------------------
p <- 0.5
b <- 1
rho <- 0.08                 # fixed edge probability: expected degree is about n*rho
n_grid <- seq(100, 2000, by = 50)
n_show <- c(200, 500, 1000, 2000)
B <- 6000                   # number of treatment-assignment draws per displayed n

# -----------------------------
# 2. Theoretical variance curve
# -----------------------------
# For G(n, rho), D_i ~ Binomial(n-1, rho), hence
# E[D_i^2] = (n-1)rho(1-rho) + (n-1)^2 rho^2.
# Since Var_Z(T_b | A) = b^2/[n^2 p(1-p)] * sum_i D_i^2,
# this is the expected conditional variance over the random graph.
expected_sum_degree_sq <- function(n, rho) {
  n * ((n - 1) * rho * (1 - rho) + (n - 1)^2 * rho^2)
}

var_exact_mean <- b^2 * expected_sum_degree_sq(n_grid, rho) /
  (n_grid^2 * p * (1 - p))

# Leading-order approximation when n*rho is large:
# Var(T_b) ~ b^2 * n * rho^2 / [p(1-p)].
var_leading <- b^2 * n_grid * rho^2 / (p * (1 - p))

# -----------------------------
# 3. Helper: simulate an ER degree sequence
# -----------------------------
er_degree_sequence <- function(n, rho) {
  degree <- integer(n)
  for (i in seq_len(n - 1L)) {
    # Edges from node i to nodes i+1,...,n
    edges <- rbinom(n - i, size = 1, prob = rho)
    degree[i] <- degree[i] + sum(edges)
    if (any(edges == 1L)) {
      neighbors <- which(edges == 1L) + i
      degree[neighbors] <- degree[neighbors] + 1L
    }
  }
  degree
}

# -----------------------------
# 4. Simulate T_b under treatment randomization
# -----------------------------
# For each n, first draw one ER graph (through its exact degree sequence).
# Then repeatedly randomize treatments, keeping that graph fixed.
# This displays the assignment variability conditional on a realized graph.
T_draws <- vector("list", length(n_show))
names(T_draws) <- as.character(n_show)
conditional_vars <- numeric(length(n_show))

for (k in seq_along(n_show)) {
  n <- n_show[k]
  deg <- er_degree_sequence(n, rho)

  # Exact conditional variance for this realized graph.
  conditional_vars[k] <- b^2 * sum(deg^2) / (n^2 * p * (1 - p))

  # Rows are independent treatment assignments; columns are individuals.
  Z <- matrix(rbinom(B * n, size = 1, prob = p), nrow = B, ncol = n)
  centered_weighted_sum <- as.vector(Z %*% deg) - p * sum(deg)
  T_draws[[k]] <- b * centered_weighted_sum / (n * p * (1 - p))
}

# -----------------------------
# 5. Produce the figure
# -----------------------------
out_file <- file.path(figure_dir, "indirect_variance_explosion.png")
png(filename = out_file, width = 2200, height = 1050, res = 180)
par(mfrow = c(1, 2), mar = c(5.0, 5.2, 3.8, 1.4),
    mgp = c(3.0, 0.9, 0), las = 1, family = "serif")

# Panel A: variance growth
plot(n_grid, var_exact_mean, type = "l", lwd = 3,
     log = "xy", xlab = "Population size n",
     ylab = expression(Var[T[b]]),
     main = "(a) The variance increases with network size",
     cex.main = 1.05, cex.lab = 1.15, cex.axis = 1.0)
lines(n_grid, var_leading, lwd = 2, lty = 2)
legend("topleft",
       legend = c("Expected conditional variance (ER graph)",
                  expression(b^2*n*rho^2/(p*(1-p)))),
       lty = c(1, 2), lwd = c(3, 2), bty = "n", cex = 0.9)
mtext(paste0("rho = ", rho, ", p = ", p, ", b = ", b),
      side = 3, line = 0.5, cex = 0.9)

# Panel B: empirical distributions for increasing n
line_types <- seq_along(n_show)
line_widths <- c(1.8, 2.1, 2.4, 2.8)
# Draw densities and determine a common vertical scale.
densities <- lapply(T_draws, density, n = 512)
y_max <- max(vapply(densities, function(d) max(d$y), numeric(1)))
plot(NULL, xlim = range(unlist(T_draws)), ylim = c(0, y_max * 1.08),
     xlab = expression(T[b]), ylab = "Density",
     main = "(b) Assignment distributions become more dispersed",
     cex.main = 1.05, cex.lab = 1.15, cex.axis = 1.0)
for (k in seq_along(densities)) {
  lines(densities[[k]], lty = line_types[k], lwd = line_widths[k])
}
legend("topright",
       legend = paste0("n = ", n_show),
       lty = line_types, lwd = line_widths,
       bty = "n", cex = 0.9)

mtext(paste0("Fixed edge probability rho = ", rho,
               "; expected degree grows approximately as n*rho"),
      side = 3, line = 0.5, cex = 0.9)

dev.off()

# -----------------------------
# 6. Print a compact numerical check
# -----------------------------
results <- data.frame(
  n = n_show,
  expected_degree = (n_show - 1) * rho,
  theoretical_mean_variance = b^2 * expected_sum_degree_sq(n_show, rho) /
    (n_show^2 * p * (1 - p)),
  realized_graph_variance = conditional_vars,
  simulated_variance = vapply(T_draws, var, numeric(1))
)
print(results, row.names = FALSE)
cat("\nFigure saved to:", normalizePath(out_file), "\n")
