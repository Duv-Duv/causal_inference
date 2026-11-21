#!/usr/bin/env Rscript
# Reproducible simulation suite for the Results section.
# Dependencies: base R, ggplot2. Optional: MASS is not required.

if (!"--only-variance" %in% commandArgs(trailingOnly = TRUE)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2: install.packages(\"ggplot2\")")
  suppressPackageStartupMessages(library(ggplot2))
}
args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
only_variance <- "--only-variance" %in% args
with_variance <- "--variance" %in% args || only_variance
root_arg <- args[!grepl("^--", args)]
# Output root: optional positional argument; otherwise current working directory.
# No pre-existing project folders are required.
root <- if (length(root_arg)) normalizePath(root_arg[1], mustWork = FALSE) else getwd()
dir.create(root, recursive = TRUE, showWarnings = FALSE)
resdir <- file.path(root, "results")
figdir <- file.path(root, "figures")
dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261125)
B <- if (quick) 30 else 180

save_csv <- function(x, name) write.csv(x, file.path(resdir, name), row.names = FALSE)
style <- theme_bw(base_size = 11, base_family = "serif") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey93"), legend.position = "bottom",
        legend.title = element_blank(), plot.title = element_text(face = "bold", size = 12),
        axis.text = element_text(colour = "grey20"))
export <- function(p, name, w = 7.2, h = 4.2) {
  ggsave(file.path(figdir, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  ggsave(file.path(figdir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
rmse <- function(x, truth) sqrt(mean((x - truth)^2, na.rm = TRUE))
summarise_mc <- function(d, keys) {
  spl <- split(d, interaction(d[keys], drop = TRUE, lex.order = TRUE))
  out <- lapply(spl, function(z) {
    truth <- z$truth[1]
    if (!is.finite(truth)) stop("Non-finite truth in Monte Carlo summary")
    estimates <- z$estimate[is.finite(z$estimate)]
    data.frame(z[1, keys, drop = FALSE], truth = truth,
               mean = mean(estimates),
               bias = mean(estimates - truth),
               sd = sd(estimates),
               rmse = sqrt(mean((estimates - truth)^2)),
               reps = length(estimates))
  })
  do.call(rbind, out)
}

# ---------- graph and causal helpers ----------
er_graph <- function(n, rho) {
  U <- matrix(runif(n*n), n, n); A <- (U < rho) * 1L
  A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0L; A
}
sbm_graph <- function(n, pin, pout) {
  block <- rep(1:2, length.out = n)
  P <- outer(block, block, function(a,b) ifelse(a == b, pin, pout))
  U <- matrix(runif(n*n), n, n); A <- (U < P) * 1L
  A[lower.tri(A)] <- t(A)[lower.tri(A)]; diag(A) <- 0L
  list(A = A, block = block)
}
exposure <- function(A, Z) {
  N <- rowSums(A); M <- as.vector(A %*% Z); X <- ifelse(N > 0, M/N, 0)
  list(N=N, M=M, X=X)
}
linear_outcome <- function(A, Z, tau=1, beta=1, gamma=1, baseline=NULL, sigma=.5) {
  e <- exposure(A,Z); n <- length(Z)
  if (is.null(baseline)) baseline <- rnorm(n, 0, sigma)
  Y <- baseline + tau*Z + beta*e$X + gamma*Z*e$X
  list(Y=Y, N=e$N, M=e$M, X=e$X)
}
finite_targets <- function(A, p=.5, tau=1, beta=1, gamma=1) {
  q <- mean(rowSums(A)>0)
  c(direct=tau+gamma*p*q, indirect=(beta+gamma*p)*q,
    total=tau+beta*q+2*gamma*p*q)
}
direct_ht <- function(Y,Z,p) mean(Y*(Z/p-(1-Z)/(1-p)))
direct_hajek <- function(Y,Z) mean(Y[Z==1]) - mean(Y[Z==0])
indirect_unbiased <- function(Y,A,Z,p) {
  e <- exposure(A,Z); g <- (e$M-p*e$N)/(p*(1-p)); mean(Y*g)
}
pc_indirect <- function(Y,A,Z,p,rank) {
  e <- exposure(A,Z); g <- (e$M-p*e$N)/(p*(1-p))
  ee <- eigen(A, symmetric=TRUE)
  ord <- order(abs(ee$values), decreasing=TRUE)[seq_len(min(rank,nrow(A)))]
  V <- ee$vectors[,ord,drop=FALSE]
  gb <- g - V %*% as.vector(crossprod(V,g))
  mean(Y*gb)
}

run_experiments <- function() {
# ============================================================
# Experiment 1. Exact exposure probabilities and positivity
# Application/reproduction: design-based exposure weighting.
# Extension: exact probabilities + density/common-support stress test.
# Complete randomization with exactly K treated units.
# Exposure labels: 00 untreated/no treated neighbor; 01 untreated/treated neighbor;
#                  10 treated/no treated neighbor; 11 treated/treated neighbor.
# ============================================================
set.seed(20261125) # Experiment 1 only; reproducible independently
choose_ratio <- function(a,b,c,d) { # choose(a,b)/choose(c,d), stable on log scale
  if (b < 0 || b > a || d < 0 || d > c) return(0)
  exp(lchoose(a,b)-lchoose(c,d))
}
exact_exposure_probs <- function(N, n, K) {
  # Own treatment is part of complete randomization; neighbors exclude self.
  p10 <- (K/n) * vapply(N, function(d) choose_ratio(n-1-d, K-1, n-1, K-1), numeric(1))
  p11 <- K/n - p10
  p00 <- ((n-K)/n) * vapply(N, function(d) choose_ratio(n-1-d, K, n-1, K), numeric(1))
  p01 <- (n-K)/n - p00
  data.frame(p00=p00,p01=p01,p10=p10,p11=p11)
}
# Experiment 1 settings (change here without affecting other experiments).
exp1_n <- 1000L
exp1_K <- 100L
exp1_densities <- c(.01, .02, .04, .06, .08, .10, .12, .14, .17, .20)
exp1_raw <- vector("list", length(exp1_densities) * B)
exp1_prob <- vector("list", length(exp1_densities))
idx <- 1L
for (j in seq_along(exp1_densities)) {
  density <- exp1_densities[j]
  for (b in seq_len(B)) {
    A <- er_graph(exp1_n, density)
    N <- rowSums(A)
    probs <- exact_exposure_probs(N, exp1_n, exp1_K)
    if (b == 1L) {
      exp1_prob[[j]] <- data.frame(density=density, degree=N, probs)
    }
    # Baseline is degree-associated; indirect exposure effect is exactly one.
    Nsd <- sd(N)
    y0 <- if (Nsd > 0) (N - mean(N)) / Nsd else rep(0, exp1_n)
    y0 <- y0 + rnorm(exp1_n, 0, .5)
    eff <- c(`00`=0, `01`=1, `10`=.5, `11`=1.5)
    Z <- integer(exp1_n)
    Z[sample.int(exp1_n, exp1_K)] <- 1L
    M <- as.vector(A %*% Z)
    E <- ifelse(Z == 0 & M == 0, "00",
                ifelse(Z == 0, "01", ifelse(M == 0, "10", "11")))
    Y <- y0 + unname(eff[E])
    # Numerical common support: excludes effectively unobservable 00 events.
    common <- probs$p01 > 1e-12 & probs$p00 > 1e-12
    if (!any(common)) stop("No common-support nodes in Experiment 1")
    # HT is defined even when an exposure group is absent in a particular draw.
    # Restrict before dividing to avoid 0/0 for nodes outside common support.
    ht <- mean(Y[common] * (as.numeric(E[common] == "01") / probs$p01[common] -
                            as.numeric(E[common] == "00") / probs$p00[common]))
    has00 <- any(E == "00" & common)
    has01 <- any(E == "01" & common)
    naive <- if (has00 && has01) {
      mean(Y[E == "01" & common]) - mean(Y[E == "00" & common])
    } else NA_real_
    exp1_raw[[idx]] <- data.frame(density=density, rep=b, truth=1,
                                   ht=ht, unweighted=naive,
                                   has00=has00, has01=has01,
                                   common_n=sum(common))
    idx <- idx + 1L
  }
}
exp1_draws <- do.call(rbind, exp1_raw)
save_csv(exp1_draws, "positivity_draws.csv")
exp1 <- rbind(
  data.frame(density=exp1_draws$density, estimator="HT exact",
             estimate=exp1_draws$ht, truth=exp1_draws$truth),
  data.frame(density=exp1_draws$density, estimator="Unweighted",
             estimate=exp1_draws$unweighted, truth=exp1_draws$truth)
)
exp1s <- summarise_mc(exp1, c("density", "estimator"))
save_csv(exp1s, "positivity_summary.csv")
probd <- do.call(rbind, exp1_prob)
save_csv(probd, "positivity_exact_probabilities.csv")

# Figure 11: one pair of exact probability curves; rho affects the observed
# degree distribution, not P(E=e | degree) at fixed n and K.
# Limit the x-axis to the degree range observed in the simulated graphs.
degree_max <- max(probd$degree)
degree_grid <- data.frame(degree=0:degree_max)
degree_grid <- cbind(degree_grid,
                     exact_exposure_probs(degree_grid$degree, exp1_n, exp1_K))
curve_data <- rbind(
  data.frame(degree=degree_grid$degree, probability=degree_grid$p01,
             exposure="01: untreated, treated neighbor"),
  data.frame(degree=degree_grid$degree, probability=degree_grid$p10,
             exposure="10: treated, no treated neighbor")
)
curve_data$exposure <- factor(curve_data$exposure,
  levels=c("01: untreated, treated neighbor", "10: treated, no treated neighbor"))
p <- ggplot(curve_data, aes(degree, probability, linetype=exposure)) +
  geom_line(linewidth=.85, colour="grey15") +
  scale_linetype_manual(values=c("solid", "longdash")) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, by=.2)) +
  labs(x="Node degree", y="Exact exposure probability",
       title="Exact exposure probabilities under complete randomization",
       subtitle=sprintf("n = %d, K = %d (10%% treated)", exp1_n, exp1_K)) +
  style
export(p, "exact_exposure_probabilities", 7.2, 4.2)

# Figure 12: retain the original two-estimator RMSE line plot. Missing
# unweighted estimates are NOT imputed; flag unreliable density settings.
validity <- aggregate(cbind(has00, has01) ~ density, exp1_draws, mean)
names(validity)[2:3] <- c("fraction_with_00", "fraction_with_01")
validity$missing00_rate <- 1 - validity$fraction_with_00
save_csv(validity, "positivity_diagnostic.csv")
exp1s <- merge(exp1s, validity, by="density", sort=FALSE)
exp1s$reliable <- exp1s$missing00_rate <= .10 &
  (exp1s$estimator != "Unweighted" | exp1s$reps >= .8 * B)
# Dashed extension for unreliable points; RMSE there is descriptive only.
p <- ggplot(exp1s, aes(density, rmse, linetype=estimator, shape=estimator,
                       group=estimator)) +
  geom_line(colour="grey25", linewidth=.65) +
  geom_point(aes(alpha=reliable), colour="grey15", size=2.3) +
  scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=.3), guide="none") +
  scale_x_continuous(breaks=exp1_densities) +
  labs(x="Edge probability", y="RMSE",
       title="Exact propensity weighting and positivity",
       subtitle="Faded points: >10% of draws lack exposure 00 (or insufficient unweighted estimates)") +
  style + theme(axis.text.x=element_text(angle=45, hjust=1),
                plot.subtitle=element_text(size=9))
export(p, "exact_positivity_rmse", 8, 4.4)
message("Experiment 1: n=", exp1_n, ", K=", exp1_K,
        "; missing-00 rates: ",
        paste(sprintf("%.2f: %.1f%%", validity$density,
                      100*validity$missing00_rate), collapse=", "))

# ============================================================
# Experiment 2. Fixed-density direct and indirect estimation
# Reproduction/application of Li--Wager mechanisms.
# ============================================================
raw2 <- list(); idx <- 1; pz <- .5
for (n in c(200,400,800)) for (rho in c(.01,.03,.06,.12)) for (b in seq_len(B)) {
  A <- er_graph(n,rho); Z <- rbinom(n,1,pz); o <- linear_outcome(A,Z); tr <- finite_targets(A,pz)
  vals <- c("Direct HT"=direct_ht(o$Y,Z,pz),"Direct Hajek"=direct_hajek(o$Y,Z),"Indirect unbiased"=indirect_unbiased(o$Y,A,Z,pz))
  truths <- c("Direct HT"=unname(tr["direct"]),
              "Direct Hajek"=unname(tr["direct"]),
              "Indirect unbiased"=unname(tr["indirect"]))
  for (m in names(vals)) {raw2[[idx]]<-data.frame(n,rho,b,estimator=m,estimate=vals[m],truth=truths[m]);idx<-idx+1}
}
r2 <- do.call(rbind,raw2); s2 <- summarise_mc(r2,c("n","rho","estimator")); save_csv(s2,"density_estimators_summary.csv")
direct <- subset(s2,estimator!="Indirect unbiased"); indirect <- subset(s2,estimator=="Indirect unbiased")
p <- ggplot(direct,aes(rho,rmse,linetype=estimator,shape=estimator,group=estimator))+geom_line()+geom_point()+facet_wrap(~n,labeller=label_both)+labs(x=expression(rho),y="RMSE",title="Direct-effect estimation under fixed density")+style; export(p,"direct_density_rmse",8,4)
p <- ggplot(indirect,aes(rho,rmse,group=1))+geom_line()+geom_point()+facet_wrap(~n,labeller=label_both)+labs(x=expression(rho),y="RMSE",title="Natural indirect estimator under fixed density")+style; export(p,"indirect_unbiased_rmse",8,4)

# ============================================================
# Experiment 3. Spectral balancing and rank sensitivity
# Reproduction/application: PC balancing on low-rank network structure.
# Extension: deliberate under/over-specification of working rank.
# ============================================================
raw3 <- list(); idx <- 1
for (n in c(300,600)) for (rho in c(.03,.06,.10)) for (b in seq_len(B)) {
  g <- sbm_graph(n, pin=min(.95,1.5*rho), pout=.5*rho); A<-g$A; Z<-rbinom(n,1,.5)
  baseline <- ifelse(g$block==1,1,-1) + rnorm(n,0,.35)
  o<-linear_outcome(A,Z,baseline=baseline); tr<-finite_targets(A,.5)["indirect"]
  e3 <- exposure(A,Z); gamma3 <- (e3$M-.5*e3$N)/(.5*.5)
  eig3 <- eigen(A, symmetric=TRUE)
  ord3 <- order(abs(eig3$values), decreasing=TRUE)[1:4]
  V3 <- eig3$vectors[,ord3,drop=FALSE]
  pc3 <- function(k) { Vk <- V3[,seq_len(k),drop=FALSE]; mean(o$Y * (gamma3 - as.vector(Vk %*% crossprod(Vk,gamma3)))) }
  vals<-c("Unbiased"=mean(o$Y*gamma3),"PC rank 1"=pc3(1),"PC rank 2"=pc3(2),"PC rank 4"=pc3(4))
  for(m in names(vals)){raw3[[idx]]<-data.frame(n,rho,b,estimator=m,estimate=unname(vals[m]),truth=unname(tr));idx<-idx+1}
}
r3<-do.call(rbind,raw3);s3<-summarise_mc(r3,c("n","rho","estimator"));save_csv(s3,"pc_balancing_summary.csv")
p<-ggplot(s3,aes(rho,rmse,linetype=estimator,shape=estimator,group=estimator))+geom_line()+geom_point()+facet_wrap(~n,labeller=label_both)+labs(x=expression(rho),y="RMSE",title="PC balancing and working-rank sensitivity")+style;export(p,"pc_rank_rmse",8,4)

# ============================================================
# Experiment 4. Sparse ER sequences, large n (100, 1000, 10000)
# Requires igraph, Matrix, RSpectra for sparse graph / rank-one eigenpair.
# Only this experiment uses sparse matrices; other experiments unchanged.
# ============================================================
if (!requireNamespace("igraph", quietly=TRUE) ||
    !requireNamespace("Matrix", quietly=TRUE) ||
    !requireNamespace("RSpectra", quietly=TRUE)) {
  stop("Experiment 4 requires packages: igraph, Matrix, RSpectra")
}
set.seed(20261125)
B4 <- if (quick) 12L else 180L
raw4 <- list(); idx <- 1L
for (aexp in c(.25,.40,.65)) for (n in c(100,1000,10000)) for (b in seq_len(B4)) {
  rho <- .55*n^(-aexp)
  graph <- igraph::sample_gnp(n, rho, directed=FALSE, loops=FALSE)
  A <- igraph::as_adjacency_matrix(graph, sparse=TRUE)
  N <- as.numeric(Matrix::rowSums(A))
  Z <- rbinom(n, 1, .5)
  M <- as.numeric(A %*% Z)
  X <- ifelse(N>0, M/N, 0)
  Y <- rnorm(n,0,.5) + Z + X + Z*X
  q <- mean(N>0)
  direct_truth <- 1 + .5*q
  indirect_truth <- 1.5*q
  gamma <- (M-.5*N)/(.5*.5)
  # Eigenvector with largest eigenvalue magnitude, matching pc_indirect().
  eig <- RSpectra::eigs_sym(A, k=1, which="LM", opts=list(tol=1e-5))
  v <- as.numeric(eig$vectors[,1])
  vals <- c("Direct HT"=mean(Y*(Z/.5-(1-Z)/.5)),
            "Direct Hajek"=mean(Y[Z==1])-mean(Y[Z==0]),
            "Indirect unbiased"=mean(Y*gamma),
            "PC rank 1"=mean(Y*(gamma-v*sum(v*gamma))))
  truths <- c("Direct HT"=direct_truth, "Direct Hajek"=direct_truth,
              "Indirect unbiased"=indirect_truth, "PC rank 1"=indirect_truth)
  for (m in names(vals)) {
    raw4[[idx]] <- data.frame(exponent=aexp,n=n,rho=rho,b=b,
        estimator=m,estimate=unname(vals[m]),truth=unname(truths[m]))
    idx <- idx+1L
  }
}
r4 <- do.call(rbind,raw4)
s4 <- summarise_mc(r4,c("exponent","n","rho","estimator"))
save_csv(s4,"scaling_summary.csv")
# Original Experiment 4 plot code retained verbatim except x-axis breaks.
p<-ggplot(s4,aes(n,rmse,linetype=estimator,shape=estimator,group=estimator))+geom_line()+geom_point(size=1.7)+facet_wrap(~exponent,labeller=label_both)+scale_x_log10(breaks=c(100,1000,10000))+scale_y_log10()+labs(x="Number of nodes (log scale)",y="RMSE (log scale)",title="Sparse-sequence diagnostic")+style;export(p,"sparsity_scaling",8,4)

# ============================================================
# Experiment 5. Network measurement error
# Extension: distinguish noisy spectral basis from noisy exposure weights.
# Edge perturbation uses deletion/addition rates scaled to graph density so that
# sparse graphs are not overwhelmed by false positives from a symmetric flip.
# ============================================================
perturb_graph <- function(A,delta) {
  n<-nrow(A); dens<-mean(A[upper.tri(A)]); del<-delta; add<-delta*dens/max(1-dens,1e-8)
  B<-A; U<-matrix(runif(n*n),n,n); upper<-upper.tri(A)
  drop<-upper & A==1 & U<del; addm<-upper & A==0 & U<add
  B[drop]<-0; B[addm]<-1
  B[lower.tri(B)] <- t(B)[lower.tri(B)]; diag(B)<-0; B
}
pc_with_basis <- function(Y,Aweights,Abasis,Z,p,rank=2) {
  e<-exposure(Aweights,Z);g<-(e$M-p*e$N)/(p*(1-p));ee<-eigen(Abasis,symmetric=TRUE);ord<-order(abs(ee$values),decreasing=TRUE)[seq_len(rank)];V<-ee$vectors[,ord,drop=FALSE];mean(Y*(g-V%*%as.vector(crossprod(V,g))))
}
raw5<-list();idx<-1
for(delta in c(0,.01,.03,.08)) for(b in seq_len(B)) {
  n<-500;g<-sbm_graph(n,.08,.02);A<-g$A;Z<-rbinom(n,1,.5);baseline<-ifelse(g$block==1,1,-1)+rnorm(n,0,.35);o<-linear_outcome(A,Z,baseline=baseline);truth<-finite_targets(A,.5)["indirect"];Anoisy<-perturb_graph(A,delta)
  vals<-c("True graph"=pc_with_basis(o$Y,A,A,Z,.5,2),"Noisy basis only"=pc_with_basis(o$Y,A,Anoisy,Z,.5,2),"Noisy basis + weights"=pc_with_basis(o$Y,Anoisy,Anoisy,Z,.5,2))
  for(m in names(vals)){raw5[[idx]]<-data.frame(error_rate=delta,b,estimator=m,estimate=vals[m],truth=truth);idx<-idx+1}
}
r5<-do.call(rbind,raw5);s5<-summarise_mc(r5,c("error_rate","estimator"));save_csv(s5,"network_error_summary.csv")
p<-ggplot(s5,aes(error_rate,rmse,linetype=estimator,shape=estimator,group=estimator))+geom_line()+geom_point()+labs(x="Edge-error severity",y="RMSE",title="Sensitivity to network measurement error")+style;export(p,"network_measurement_error",8,4)

message("Completed ", B, " Monte Carlo replications per configuration.")
message("Results: ", resdir)
message("Figures: ", figdir)

}

if (!only_variance) run_experiments()

# Standalone illustration: run with --variance or --only-variance.
run_variance_illustration <- function() {
  set.seed(2026)
  variance_figure_dir <- file.path(root, "figures", "academic")
  dir.create(variance_figure_dir, recursive = TRUE, showWarnings = FALSE)
# -----------------------------
# 1. Parameters
# -----------------------------
p <- 0.5
b <- 1
rho <- 0.08                 # fixed edge probability: expected degree is about n*rho
n_grid <- seq(100, 2000, by = 50)
n_show <- c(200, 500, 1000, 2000)
B <- if (quick) 300L else 6000L                   # number of treatment-assignment draws per displayed n

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
out_file <- file.path(variance_figure_dir, "indirect_variance_explosion.png")
draw_variance_figure <- function() {

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


}
png(filename = out_file, width = 2200, height = 1050, res = 180)
draw_variance_figure()
dev.off()
pdf(file.path(variance_figure_dir, "indirect_variance_explosion.pdf"), width = 12.2, height = 5.8, family = "serif", useDingbats = FALSE)
draw_variance_figure()
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

}
if (with_variance) run_variance_illustration()
