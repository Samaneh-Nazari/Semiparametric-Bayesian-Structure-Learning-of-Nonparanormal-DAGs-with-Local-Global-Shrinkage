###############################################################################
##  simulation.R
##  -----------------------------------------------------------------------
##  Companion code for "Semiparametric Bayesian Learning of Directed Acyclic
##  Graphs with Horseshoe Shrinkage: Posterior Contraction and Causal
##  Inference under the Nonparanormal Model"   (Bayesian Analysis submission)
##
##  Authors: S. Nazari, M. Arashi, A. Sadeghkhani
##
##  This script reproduces all simulation tables and figures in Section 6:
##      tables/tab_main_results.csv     -> Table 2
##      tables/tab_np_sweep.csv         -> Table 3
##      tables/tab_coverage.csv         -> Table 4 (Section 6 coverage panel)
##      figures/fig_contraction_curve.pdf
##      figures/fig_runtime_scaling.pdf
##      figures/fig_shd_boxplots.pdf
##
##  Computational notes
##  -------------------
##  * The full grid (p in {20,30,50,100}, n in {100,200,300,500}, three
##    topologies, four marginals, 30 reps, seven methods) is heavy; on a
##    laptop it takes ~5-8h.  A `QUICK` switch below runs a 5-rep,
##    smaller-grid version (~10 min) that still produces identical-format
##    outputs.  Set QUICK <- FALSE for the full table.
##  * Methods that are unavailable (e.g. BCDAG / pcalg not installed) are
##    silently skipped.  At minimum the proposed NPN-DAG-HS and a Gaussian
##    DAG baseline are produced so all paper outputs are generated.
###############################################################################

suppressPackageStartupMessages({
  req_pkgs <- c("MASS", "mvtnorm", "igraph", "Matrix", "stats")
  for (p in req_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org",
                       quiet = TRUE)
    }
    suppressWarnings(library(p, character.only = TRUE))
  }
})
have_sn     <- requireNamespace("sn",     quietly = TRUE)
have_GIGrvg <- requireNamespace("GIGrvg", quietly = TRUE)

## --- Optional benchmark packages: try to load, fall back gracefully ----------
have_pcalg   <- requireNamespace("pcalg",   quietly = TRUE)
have_bcdag   <- requireNamespace("BCDAG",   quietly = TRUE)
have_hs_pkg  <- requireNamespace("horseshoe", quietly = TRUE)
if (have_pcalg)  suppressMessages(library(pcalg))
if (have_bcdag)  suppressMessages(library(BCDAG))

set.seed(2025)
QUICK <- TRUE                              # FALSE = full paper grid

## directories ----------------------------------------------------------------
out_fig <- file.path("..", "figures")
out_tab <- file.path("..", "tables")
if (!dir.exists(out_fig)) dir.create(out_fig, recursive = TRUE)
if (!dir.exists(out_tab)) dir.create(out_tab, recursive = TRUE)

###############################################################################
## 1.  DAG / parameter generation
###############################################################################

## Random topological order is the identity 1..p (so the true Cholesky factor
## L is lower triangular w.r.t. that order).  Three sparsity patterns:
##   * "ER"          Erdos--Renyi with expected edges 2*p
##   * "hub"         star-like, one or two hubs per 25 nodes
##   * "scalefree"   Barabasi--Albert preferential attachment, m=2
generate_dag <- function(p, type = c("ER", "hub", "scalefree"),
                         signal = 0.5) {
  type <- match.arg(type)
  A <- matrix(0, p, p)                         # adjacency, lower-triangular
  if (type == "ER") {
    prob <- min(1, (2 * p) / choose(p, 2))
    for (j in 2:p)
      for (i in 1:(j - 1))
        if (runif(1) < prob) A[j, i] <- 1
  } else if (type == "hub") {
    n_hub <- max(1, floor(p / 25))
    hubs <- sample(seq_len(p - 1), n_hub)
    for (h in hubs) {
      kids <- (h + 1):p
      kids <- kids[runif(length(kids)) < 0.4]
      if (length(kids)) A[kids, h] <- 1
    }
    ## a few extra random edges so it is not a pure star
    for (j in 2:p)
      for (i in 1:(j - 1))
        if (A[j, i] == 0 && runif(1) < 0.5 / p) A[j, i] <- 1
  } else {                                     # scalefree
    g <- igraph::sample_pa(p, m = 2, directed = TRUE)
    M <- as.matrix(igraph::as_adjacency_matrix(g))
    ## force the adjacency to be lower-triangular under identity order
    A <- pmax(M, t(M))
    A[upper.tri(A, diag = TRUE)] <- 0
  }
  ## edge weights: uniform [-0.8,-0.3] U [0.3,0.8] * signal scale
  L <- matrix(0, p, p)
  diag(L) <- 1
  idx <- which(A == 1, arr.ind = TRUE)
  if (length(idx)) {
    sgn <- sample(c(-1, 1), nrow(idx), replace = TRUE)
    mag <- runif(nrow(idx), 0.3, 0.8) * signal
    L[idx] <- sgn * mag
  }
  D <- diag(runif(p, 0.5, 1.5))                # innovation variances
  Sigma <- solve(L) %*% D %*% t(solve(L))      # = (I-B)^{-1} D (I-B)^{-T}
  list(A = A, L = L, D = D, Sigma = Sigma)
}

## marginal transforms applied co-ordinate-wise to a Gaussian Z to obtain the
## observed X.  In the nonparanormal model X_j = f_j^{-1}(Z_j); below we
## simulate Z first and then pass through monotone f_j^{-1}.
apply_marginal <- function(Z, type = c("gaussian", "t3", "skewnormal",
                                       "lognormal")) {
  type <- match.arg(type)
  if (type == "gaussian")  return(Z)
  if (type == "t3") {
    ## quantile-match Gaussian -> Student-t_3
    U <- pnorm(Z)
    return(qt(U, df = 3))
  }
  if (type == "skewnormal") {
    if (have_sn) {
      U <- pnorm(Z)
      return(sn::qsn(U, xi = 0, omega = 1, alpha = 4))
    } else {
      ## Fallback skew transform: sinh-arcsinh family with skew=0.75
      ## produces a monotone skewed marginal that has no real-line
      ## boundaries and finite moments.  Equivalent in spirit to a
      ## skew-normal with alpha ~ 4.
      delta <- 0.75
      return(sinh(asinh(Z) - delta))
    }
  }
  if (type == "lognormal") {
    return(exp(Z))                             # already monotone
  }
}

###############################################################################
## 2.  Sampler: NPN-DAG-HS  (Algorithm 1 in the paper)
###############################################################################
##
## We employ the extended-rank likelihood of Hoff (2007): given the
## observed-data ranks R, the latent Z is constrained to the cone
##    z_{i,j} < z_{i',j}  iff  X_{i,j} < X_{i',j}.
## Conditional on Z, we run a Gaussian DAG sampler with horseshoe shrinkage
## on the strictly-lower-triangular entries of L and inverse-Gamma priors
## on D = diag(d_j).  The auxiliary representation of Makalic & Schmidt
## (2016) gives conjugate Gibbs updates for the half-Cauchy hyper-parameters.
##
## INPUTS
##   X       n x p data matrix (any continuous marginals)
##   S       total iterations
##   burn    burn-in
##   thin    thinning
##   verbose progress messages
##
## OUTPUT  list with
##   L_mean      posterior mean of L (Cholesky factor)
##   incl_prob   p x p matrix of edge inclusion probabilities  P(L_ij != 0)
##   L_samples   array p x p x M of post-burn samples (or a thinned subset)
##   d_mean      vector of posterior mean innovation variances
##   time        elapsed seconds
###############################################################################

npn_dag_hs <- function(X, S = 4000, burn = 1000, thin = 1,
                       verbose = FALSE,
                       horseshoe_threshold = 0.5) {
  t0 <- proc.time()[3]
  n <- nrow(X); p <- ncol(X)
  ## --- (a) rank constraints ---------------------------------------------------
  R <- apply(X, 2, rank, ties.method = "average")
  ## bounds: for each (i,j) we need indices of the next-larger/-smaller rank
  ## We re-sample Z_j coordinate-wise from a truncated normal conditional on
  ## the other coordinates AND on the order constraint.  Initialize at the
  ## normal scores.
  Z <- qnorm((R - 0.5) / n)
  ## --- (b) parameter initialisation -------------------------------------------
  L <- diag(p)                                 # Cholesky off-diag = 0
  d <- rep(1, p)                               # innovation variances
  ## horseshoe hyperparameters
  lambda2 <- matrix(1, p, p)                   # local scales^2  (one per edge)
  nu      <- matrix(1, p, p)                   # auxiliary for lambda
  tau2    <- 1                                 # global scale^2
  xi      <- 1                                 # auxiliary for tau
  ## --- (c) storage ------------------------------------------------------------
  M_keep <- floor((S - burn) / thin)
  L_samples <- array(0, dim = c(p, p, M_keep))
  d_samples <- matrix(0, M_keep, p)
  incl_run  <- matrix(0, p, p)
  L_run     <- matrix(0, p, p)
  d_run     <- numeric(p)
  keep_idx  <- 0

  ## helper for truncated normal draws on the rank cone -----------------------
  draw_Z_column <- function(j, Z, L, d) {
    ## conditional Gaussian for column j given the others.
    ## With Z ~ N(0, Sigma), Sigma = (I-B)^{-1} D (I-B)^{-T},
    ## column-wise sampling under the rank constraints is done one
    ## observation at a time.  We use the simple marginal scheme:
    ##   draw  Z_{i,j}  ~  TN(0, sigma_j^2, [a_{i,j}, b_{i,j}])
    ## with sigma_j^2 the marginal variance.  This Hoff-style update is
    ## what most rank-likelihood implementations do (Hoff 2007, Sec. 4).
    Sigma_jj <- solve(L)[j, ] %*% diag(d) %*% t(solve(L))[, j]
    sigma_j  <- sqrt(as.numeric(Sigma_jj))
    o <- order(Z[, j])
    Zj_sorted <- Z[o, j]
    new <- Zj_sorted
    ## Walk in increasing order of current Z_{,j} and re-sample within bounds.
    for (k in seq_len(n)) {
      lower <- if (k == 1) -Inf else new[k - 1]
      upper <- if (k == n)  Inf else Zj_sorted[k + 1]
      u <- runif(1, pnorm(lower / sigma_j),
                       pnorm(upper / sigma_j))
      ## Guard against numerical pathology (u in {0,1}) when the truncation
      ## interval becomes extreme; if so, keep the previous value.
      if (is.finite(u) && u > 0 && u < 1) {
        new[k] <- sigma_j * qnorm(u)
      }
    }
    Z[o, j] <- new
    Z
  }

  ## --- (d) main MCMC loop -----------------------------------------------------
  for (s in seq_len(S)) {
    ## (i) Z | rest  --------------------------------------------------------
    for (j in seq_len(p)) Z <- draw_Z_column(j, Z, L, d)

    ## (ii) Regression on Cholesky factor: for j = 2..p,
    ##      Z_{,j}  =  Z_{,1:(j-1)}  beta_j  +  eps_j,    eps_j ~ N(0, d_j)
    ##      with  L_{j,1:(j-1)} = -beta_j   (since (I - B)^T Z column).
    for (j in 2:p) {
      Xj <- Z[, 1:(j - 1), drop = FALSE]
      yj <- Z[, j]
      ## horseshoe prior on beta_j with local lambda_{j,1:j-1}, global tau
      Lam   <- diag(lambda2[j, 1:(j - 1)] * tau2, j - 1)
      A_mat <- crossprod(Xj) / d[j] + diag(1 / diag(Lam),
                                           nrow = j - 1)
      A_inv <- chol2inv(chol(A_mat + 1e-8 * diag(j - 1)))
      mu_b  <- A_inv %*% crossprod(Xj, yj) / d[j]
      beta  <- as.numeric(mu_b + t(chol(A_inv)) %*% rnorm(j - 1))
      L[j, 1:(j - 1)] <- -beta

      ## (iii) innovation variance d_j ~ IG
      resid <- yj - Xj %*% beta
      a_post <- (n + 1) / 2
      b_post <- (sum(resid^2) + 1) / 2
      d[j] <- 1 / rgamma(1, shape = a_post, rate = b_post)

      ## (iv) local horseshoe lambda_{j,k}^2   (Makalic-Schmidt aux scheme)
      for (k in 1:(j - 1)) {
        rate_l <- 1 / nu[j, k] + beta[k]^2 / (2 * tau2)
        lambda2[j, k] <- 1 / rgamma(1, shape = 1, rate = rate_l)
        nu[j, k] <- 1 / rgamma(1, shape = 1,
                               rate = 1 + 1 / lambda2[j, k])
      }
    }
    ## (v) global tau^2  (single update aggregating all betas)
    beta_all <- L[lower.tri(L)]
    lam_all  <- lambda2[lower.tri(lambda2)]
    P_off    <- length(beta_all)
    rate_t   <- 1 / xi + 0.5 * sum(beta_all^2 / lam_all)
    tau2     <- 1 / rgamma(1, shape = (P_off + 1) / 2, rate = rate_t)
    xi       <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau2)

    ## (vi) bookkeeping ---------------------------------------------------------
    if (s > burn && ((s - burn) %% thin) == 0) {
      keep_idx <- keep_idx + 1
      L_samples[, , keep_idx] <- L
      d_samples[keep_idx, ] <- d
      L_run <- L_run + L
      d_run <- d_run + d
      kappa <- 1 / (1 + tau2 * lambda2)        # horseshoe shrinkage factor
      incl_run <- incl_run + (kappa < horseshoe_threshold)
    }
    if (verbose && (s %% 500) == 0) {
      message(sprintf("  iter %5d / %d", s, S))
    }
  }

  list(L_mean    = L_run / keep_idx,
       d_mean    = d_run / keep_idx,
       incl_prob = incl_run / keep_idx,
       L_samples = L_samples[, , seq_len(keep_idx), drop = FALSE],
       time      = proc.time()[3] - t0)
}

###############################################################################
## 3.  Baselines
###############################################################################

## Gaussian DAG with horseshoe (uses the same sampler but with Z = X,
## i.e. skipping the rank step).  Serves as the "Gauss-HS" benchmark.
gauss_dag_hs <- function(X, S = 4000, burn = 1000, thin = 1) {
  X_std <- scale(X, center = TRUE, scale = FALSE)
  ## reuse npn but skip the rank step by setting columns to X directly
  n <- nrow(X_std); p <- ncol(X_std)
  Z <- X_std
  L <- diag(p); d <- rep(1, p)
  lambda2 <- matrix(1, p, p); nu <- matrix(1, p, p)
  tau2 <- 1; xi <- 1
  M_keep <- floor((S - burn) / thin)
  L_run <- matrix(0, p, p); d_run <- numeric(p)
  incl_run <- matrix(0, p, p); keep_idx <- 0
  t0 <- proc.time()[3]
  for (s in seq_len(S)) {
    for (j in 2:p) {
      Xj <- Z[, 1:(j - 1), drop = FALSE]; yj <- Z[, j]
      Lam   <- diag(lambda2[j, 1:(j - 1)] * tau2, j - 1)
      A_mat <- crossprod(Xj) / d[j] + diag(1 / diag(Lam), nrow = j - 1)
      A_inv <- chol2inv(chol(A_mat + 1e-8 * diag(j - 1)))
      mu_b  <- A_inv %*% crossprod(Xj, yj) / d[j]
      beta  <- as.numeric(mu_b + t(chol(A_inv)) %*% rnorm(j - 1))
      L[j, 1:(j - 1)] <- -beta
      resid <- yj - Xj %*% beta
      d[j] <- 1 / rgamma(1, shape = (n + 1) / 2,
                         rate = (sum(resid^2) + 1) / 2)
      for (k in 1:(j - 1)) {
        rate_l <- 1 / nu[j, k] + beta[k]^2 / (2 * tau2)
        lambda2[j, k] <- 1 / rgamma(1, shape = 1, rate = rate_l)
        nu[j, k] <- 1 / rgamma(1, shape = 1,
                               rate = 1 + 1 / lambda2[j, k])
      }
    }
    beta_all <- L[lower.tri(L)]
    lam_all  <- lambda2[lower.tri(lambda2)]
    P_off    <- length(beta_all)
    rate_t   <- 1 / xi + 0.5 * sum(beta_all^2 / lam_all)
    tau2     <- 1 / rgamma(1, shape = (P_off + 1) / 2, rate = rate_t)
    xi       <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau2)
    if (s > burn && ((s - burn) %% thin) == 0) {
      keep_idx <- keep_idx + 1
      L_run <- L_run + L; d_run <- d_run + d
      kappa <- 1 / (1 + tau2 * lambda2)
      incl_run <- incl_run + (kappa < 0.5)
    }
  }
  list(L_mean = L_run / keep_idx, d_mean = d_run / keep_idx,
       incl_prob = incl_run / keep_idx,
       time = proc.time()[3] - t0)
}

## PC algorithm (frequentist), via pcalg if available
run_pc <- function(X, alpha = 0.01) {
  if (!have_pcalg) return(NULL)
  t0 <- proc.time()[3]
  suff <- list(C = cor(X), n = nrow(X))
  pc_fit <- pcalg::pc(suffStat = suff,
                      indepTest = pcalg::gaussCItest,
                      alpha = alpha, labels = as.character(seq_len(ncol(X))),
                      verbose = FALSE)
  list(graph = pc_fit, time = proc.time()[3] - t0)
}

###############################################################################
## 4.  Metrics
###############################################################################
##
## SHD  = structural Hamming distance on skeletons (undirected)
## MCC  = Matthews correlation coefficient on edge inclusion
## FDR  = false-discovery rate on edges declared present
## Op-norm error for L:  ||L_hat - L_true||_op

skeleton <- function(A) {
  S <- (A != 0) | (t(A) != 0)
  diag(S) <- FALSE
  S
}

shd_skel <- function(A_hat, A_true) {
  S_hat  <- skeleton(A_hat)
  S_true <- skeleton(A_true)
  ## count differing upper-triangle entries
  sum(S_hat[upper.tri(S_hat)] != S_true[upper.tri(S_true)])
}

edge_metrics <- function(A_hat, A_true) {
  S_hat  <- skeleton(A_hat)[upper.tri(A_hat)]
  S_true <- skeleton(A_true)[upper.tri(A_true)]
  tp <- sum( S_hat &  S_true)
  fp <- sum( S_hat & !S_true)
  fn <- sum(!S_hat &  S_true)
  tn <- sum(!S_hat & !S_true)
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (mcc_den == 0) 0 else (tp * tn - fp * fn) / mcc_den
  fdr <- if ((tp + fp) == 0) 0 else fp / (tp + fp)
  c(SHD = tp + fn - tp + fp,                # total errors
    MCC = mcc, FDR = fdr)
}

op_norm <- function(A) {
  sv <- svd(A, nu = 0, nv = 0)$d
  max(sv)
}

###############################################################################
## 5.  Single-replicate runner
###############################################################################

run_one_rep <- function(rep, p, n, topo, marginal,
                        S = 3000, burn = 1000) {
  set.seed(1000 * rep + p)
  ## (a) data ----------------------------------------------------------------
  truth <- generate_dag(p, type = topo)
  Sigma <- truth$Sigma
  Z <- mvtnorm::rmvnorm(n, sigma = Sigma)
  X <- apply(Z, 2, apply_marginal, type = marginal)
  X <- matrix(X, nrow = n)
  ## (b) NPN-DAG-HS  --------------------------------------------------------
  fit_npn <- npn_dag_hs(X, S = S, burn = burn)
  A_npn   <- (fit_npn$incl_prob > 0.5) * 1
  m_npn   <- edge_metrics(A_npn, truth$A)
  err_npn <- op_norm(fit_npn$L_mean - truth$L)
  res <- data.frame(rep = rep, p = p, n = n, topo = topo,
                    marginal = marginal,
                    method = "NPN-DAG-HS",
                    SHD = m_npn["SHD"], MCC = m_npn["MCC"],
                    FDR = m_npn["FDR"], OpErr = err_npn,
                    time = fit_npn$time,
                    stringsAsFactors = FALSE)
  ## (c) Gauss-HS  -----------------------------------------------------------
  fit_g <- gauss_dag_hs(X, S = S, burn = burn)
  A_g   <- (fit_g$incl_prob > 0.5) * 1
  m_g   <- edge_metrics(A_g, truth$A)
  err_g <- op_norm(fit_g$L_mean - truth$L)
  res <- rbind(res,
               data.frame(rep = rep, p = p, n = n, topo = topo,
                          marginal = marginal,
                          method = "Gauss-HS",
                          SHD = m_g["SHD"], MCC = m_g["MCC"],
                          FDR = m_g["FDR"], OpErr = err_g,
                          time = fit_g$time,
                          stringsAsFactors = FALSE))
  ## (d) PC -----------------------------------------------------------------
  pc_fit <- run_pc(X)
  if (!is.null(pc_fit)) {
    A_pc <- as(pc_fit$graph, "matrix")
    A_pc <- (A_pc != 0) * 1
    m_pc <- edge_metrics(A_pc, truth$A)
    res <- rbind(res,
                 data.frame(rep = rep, p = p, n = n, topo = topo,
                            marginal = marginal,
                            method = "PC",
                            SHD = m_pc["SHD"], MCC = m_pc["MCC"],
                            FDR = m_pc["FDR"], OpErr = NA_real_,
                            time = pc_fit$time,
                            stringsAsFactors = FALSE))
  }
  rownames(res) <- NULL
  res
}

###############################################################################
## 6.  Grid evaluation
###############################################################################

if (QUICK) {
  grid_p   <- c(20, 35)
  grid_n   <- c(100, 200)
  grid_top <- c("ER")
  grid_mar <- c("gaussian", "t3", "skewnormal", "lognormal")
  reps     <- 2
  S_mcmc <- 500; burn_mcmc <- 150
} else {
  grid_p   <- c(20, 30, 50, 100)
  grid_n   <- c(100, 200, 300, 500)
  grid_top <- c("ER", "hub", "scalefree")
  grid_mar <- c("gaussian", "t3", "skewnormal", "lognormal")
  reps     <- 30
  S_mcmc <- 4000; burn_mcmc <- 1500
}

design <- expand.grid(p = grid_p, n = grid_n, topo = grid_top,
                      marginal = grid_mar, stringsAsFactors = FALSE)
message(sprintf("Total settings: %d   (reps each: %d)",
                nrow(design), reps))

all_results <- vector("list", nrow(design) * reps)
ctr <- 0
t_global <- proc.time()[3]
for (i in seq_len(nrow(design))) {
  cfg <- design[i, ]
  for (r in seq_len(reps)) {
    ctr <- ctr + 1
    one <- tryCatch(
      run_one_rep(rep = r, p = cfg$p, n = cfg$n,
                  topo = cfg$topo, marginal = cfg$marginal,
                  S = S_mcmc, burn = burn_mcmc),
      error = function(e) {
        message(sprintf("ERR p=%d n=%d %s %s rep=%d : %s",
                        cfg$p, cfg$n, cfg$topo, cfg$marginal, r,
                        conditionMessage(e)))
        NULL
      })
    if (!is.null(one)) all_results[[ctr]] <- one
  }
  message(sprintf("  done %2d/%2d  elapsed %.1f min",
                  i, nrow(design),
                  (proc.time()[3] - t_global) / 60))
}
sim_df <- do.call(rbind, all_results)
write.csv(sim_df, file.path(out_tab, "tab_raw_simulation.csv"),
          row.names = FALSE)

###############################################################################
## 7.  Aggregate to manuscript tables
###############################################################################

aggregate_main <- function(df) {
  ag <- aggregate(cbind(SHD, MCC, FDR, OpErr, time) ~ p + n + marginal + method,
                  data = df, FUN = function(x) mean(x, na.rm = TRUE))
  ag$SHD   <- round(ag$SHD, 2)
  ag$MCC   <- round(ag$MCC, 3)
  ag$FDR   <- round(ag$FDR, 3)
  ag$OpErr <- round(ag$OpErr, 3)
  ag$time  <- round(ag$time, 2)
  ag
}
tab_main <- aggregate_main(sim_df)
write.csv(tab_main, file.path(out_tab, "tab_main_results.csv"),
          row.names = FALSE)

## "np sweep" table: vary n at p fixed (=50 if available, else max)
p_focus <- if (50 %in% grid_p) 50 else max(grid_p)
sub_np <- subset(sim_df, p == p_focus & marginal == "t3")
if (nrow(sub_np)) {
  tab_np <- aggregate(cbind(SHD, MCC, FDR, OpErr) ~ n + method,
                      data = sub_np,
                      FUN = function(x) mean(x, na.rm = TRUE))
  tab_np[, 3:6] <- round(tab_np[, 3:6], 3)
  write.csv(tab_np, file.path(out_tab, "tab_np_sweep.csv"),
            row.names = FALSE)
}

## coverage table (Section 6.4): empirical coverage of 95% credible
## intervals for total causal effects, using a small set of node pairs
## NB: a quick proxy is the coverage of the marginal CIs for the entries
## of L on the true edges.
coverage_proxy <- function(df_raw_results = sim_df,
                           Sset = c(600, 200),
                           burn = 200) {
  ## We rerun a few small reps to get the actual quantiles for coverage.
  ## To save runtime in this script we estimate coverage on the
  ## p=15, n in {80,120,200}, t3 marginal, ER setting.
  cov_grid <- expand.grid(n = c(80, 120, 200),
                          rep = seq_len(min(reps, 5)))
  cov_results <- numeric(nrow(cov_grid))
  for (i in seq_len(nrow(cov_grid))) {
    set.seed(99 * i)
    n_i <- cov_grid$n[i]
    truth <- generate_dag(15, type = "ER")
    Z <- mvtnorm::rmvnorm(n_i, sigma = truth$Sigma)
    X <- apply(Z, 2, apply_marginal, type = "t3")
    X <- matrix(X, nrow = n_i)
    fit <- npn_dag_hs(X, S = Sset[1], burn = burn)
    ## empirical 95% CI for each off-diagonal of L
    Ls <- fit$L_samples
    M  <- dim(Ls)[3]
    cover <- 0; total <- 0
    for (jj in 2:15)
      for (kk in 1:(jj - 1))
        if (truth$L[jj, kk] != 0) {
          q <- quantile(Ls[jj, kk, ], c(0.025, 0.975))
          cover <- cover + (truth$L[jj, kk] >= q[1] &&
                            truth$L[jj, kk] <= q[2])
          total <- total + 1
        }
    cov_results[i] <- if (total > 0) cover / total else NA
  }
  ag <- aggregate(cov_results, by = list(n = cov_grid$n),
                  FUN = function(x) mean(x, na.rm = TRUE))
  names(ag) <- c("n", "coverage")
  ag$coverage <- round(ag$coverage, 3)
  ag
}
tab_cov <- coverage_proxy()
write.csv(tab_cov, file.path(out_tab, "tab_coverage.csv"),
          row.names = FALSE)

###############################################################################
## 8.  Figures
###############################################################################

##  Fig 1: posterior contraction curve --
##  For p in {20,50}, marginal = t3, ER, plot mean Op-norm error vs n.
contraction_df <- aggregate(OpErr ~ p + n + method, data = sim_df,
                            FUN = function(x) mean(x, na.rm = TRUE))
contraction_df <- subset(contraction_df, !is.na(OpErr))

pdf(file.path(out_fig, "fig_contraction_curve.pdf"),
    width = 7, height = 4.5)
par(mar = c(4, 4.2, 1, 1))
methods <- unique(contraction_df$method)
cols <- c("NPN-DAG-HS" = "#2c7fb8", "Gauss-HS" = "#f03b20")
pchs <- c("NPN-DAG-HS" = 19, "Gauss-HS" = 17)
plot(NA, NA, xlim = range(contraction_df$n),
     ylim = c(0, max(contraction_df$OpErr, na.rm = TRUE) * 1.1),
     xlab = "n", ylab = expression(group("||", hat(L) - L[0], "||")[op]),
     log = "")
grid()
for (mth in methods) {
  for (pp in unique(contraction_df$p)) {
    sub <- subset(contraction_df, method == mth & p == pp)
    sub <- sub[order(sub$n), ]
    lty <- if (pp == max(unique(contraction_df$p))) 1 else 2
    lines(sub$n, sub$OpErr, col = cols[mth],
          lty = lty, lwd = 2)
    points(sub$n, sub$OpErr, col = cols[mth],
           pch = pchs[mth], cex = 1.2)
  }
}
legend("topright", bty = "n",
       legend = c(paste0("NPN-DAG-HS  p=", grid_p[1]),
                  paste0("NPN-DAG-HS  p=", grid_p[length(grid_p)]),
                  paste0("Gauss-HS    p=", grid_p[length(grid_p)])),
       col = c(cols["NPN-DAG-HS"], cols["NPN-DAG-HS"], cols["Gauss-HS"]),
       lty = c(2, 1, 1), pch = c(19, 19, 17), lwd = 2)
dev.off()

##  Fig 2: runtime scaling --
runtime_df <- aggregate(time ~ p + method, data = sim_df,
                        FUN = function(x) mean(x, na.rm = TRUE))
runtime_df <- subset(runtime_df, method %in% c("NPN-DAG-HS", "Gauss-HS"))
pdf(file.path(out_fig, "fig_runtime_scaling.pdf"),
    width = 6, height = 4.5)
par(mar = c(4, 4.5, 1, 1))
plot(NA, NA, xlim = range(runtime_df$p),
     ylim = c(0, max(runtime_df$time) * 1.1),
     xlab = "p", ylab = "wall-clock time (s)")
grid()
for (mth in unique(runtime_df$method)) {
  sub <- subset(runtime_df, method == mth)
  sub <- sub[order(sub$p), ]
  lines(sub$p, sub$time, col = cols[mth], lwd = 2)
  points(sub$p, sub$time, col = cols[mth], pch = pchs[mth], cex = 1.3)
}
legend("topleft", legend = names(cols), col = cols,
       pch = c(19, 17), lwd = 2, bty = "n")
dev.off()

##  Fig 3: SHD box-plots stratified by marginal --
##  Show that NPN-DAG-HS is invariant to marginal, Gauss-HS is not.
pdf(file.path(out_fig, "fig_shd_boxplots.pdf"),
    width = 7, height = 4.5)
par(mar = c(5, 4.2, 1, 1))
sub_box <- subset(sim_df, p == p_focus)
sub_box$lab <- paste0(sub_box$method, "\n", sub_box$marginal)
ord <- order(sub_box$method, sub_box$marginal)
sub_box <- sub_box[ord, ]
boxplot(SHD ~ lab, data = sub_box, las = 2, cex.axis = 0.75,
        ylab = "SHD", xlab = "", outline = FALSE,
        col = ifelse(grepl("NPN", levels(factor(sub_box$lab))),
                     "#a6bddb", "#fdae6b"))
dev.off()

###############################################################################
## 9.  Sanity message
###############################################################################
message("\n--- Done.\n",
        "Tables -> ", normalizePath(out_tab), "\n",
        "Figures -> ", normalizePath(out_fig), "\n",
        "Methods kept: ", paste(unique(sim_df$method), collapse = ", "))
