# Shared specification and implementation for Appendix B.1--B.2.
# Keep this file beside the four simulation drivers on other devices.
# No simulation is run when this file is sourced.

sim_settings <- c("1trial_binary_main", "5trial_binary_main",
                  "1trial_continuous_add", "5trial_binary_small_add")

sim_grid <- function(setting) {
  stopifnot(setting %in% sim_settings)
  single <- startsWith(setting, "1trial")
  if (single) {
    d <- expand.grid(N = c(1000L, 2000L), eta = 1:2, beta = 1:2,
                     theta = 1:2, KEEP.OUT.ATTRS = FALSE)
    d$gamma <- NA_integer_
  } else {
    d <- expand.grid(N = if (setting == "5trial_binary_main") 5000L else 2000L,
                     eta = 1:2, beta = 1:2, gamma = 1:2, theta = 1:2,
                     KEEP.OUT.ATTRS = FALSE)
  }
  d <- cbind(scenario = seq_len(nrow(d)), d)
  d[, c("scenario", "N", "eta", "beta", "gamma", "theta")]
}

# Gaussian quadrature for independent standard normals; no extra dependency.
sim_normal_rule <- function(n = 60L) {
  j <- matrix(0, n, n)
  j[cbind(seq_len(n - 1L), 2:n)] <- sqrt(seq_len(n - 1L))
  e <- eigen(j + t(j), symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}

# Integrate out the two binary variables conditional on a common normal U.
# X3 = mu + sqrt(rho)*U + sqrt(1-rho)*E3.
sim_integration_grid <- function(eta, order = 60L) {
  q <- sim_normal_rule(order)
  z <- expand.grid(u = q$x, e = q$x)
  p1 <- pnorm((sqrt(eta[4]) * z$u - qnorm(1 - eta[1])) / sqrt(1 - eta[4]))
  p2 <- pnorm((sqrt(eta[4]) * z$u - qnorm(1 - eta[2])) / sqrt(1 - eta[4]))
  w <- as.vector(outer(q$w, q$w))
  parts <- lapply(0:3, function(k) {
    x1 <- k %% 2; x2 <- k %/% 2
    data.frame(X1 = x1, X2 = x2,
      X3 = eta[3] + sqrt(eta[4]) * z$u + sqrt(1 - eta[4]) * z$e,
      w = w * if (x1) p1 else (1 - p1))
  })
  for (k in seq_along(parts)) parts[[k]]$w <- parts[[k]]$w *
    if (parts[[k]]$X2[1]) p2 else (1 - p2)
  do.call(rbind, parts)
}

sim_basis <- function(x) cbind(`(Intercept)` = 1, as.matrix(x[, c("X1", "X2", "X3")]))

sim_parameters <- function(setting, scenario, order = 60L) {
  row <- sim_grid(setting)
  stopifnot(length(scenario) == 1L, scenario %in% row$scenario)
  row <- row[row$scenario == scenario, ]
  single <- startsWith(setting, "1trial")
  small <- setting == "5trial_binary_small_add"
  eta <- list(c(.3, .3, 0, .3), c(.5, .5, 0, .5))[[row$eta]]
  beta <- if (single) {
    list(log(c(1.2, .5, .5, .5)), log(c(.5, 2, 2, 2)))[[row$beta]]
  } else {
    list(log(c(2, .5, .5, .5)), log(c(.8, 2, 2, 2)))[[row$beta]]
  }
  theta <- list(rbind(log(c(.5, 2, .5, 1.25)), log(c(.5, .5, 2, .8))),
                rbind(log(c(.5, 1.25, .8, 1.1)), log(c(.5, .8, 1.25, .9))))[[row$theta]]
  gamma <- matrix(0, if (single) 1L else 5L, 4L)
  if (!single) {
    a <- log(1.25)
    gamma[2:5, ] <- cbind(if (small) 0 else log(2),
      rbind(c(-a, a, -a), c(-a/2, a/2, -a/2),
            c(-a, a/2, a/2), c(-a/2, a, a/2)))
    if (row$gamma == 2L) gamma[, 2:4] <- -gamma[, 2:4]
  }
  integration <- sim_integration_grid(eta, order)
  b <- sim_basis(integration)
  if (small) {
    slope <- as.vector(b[, -1L] %*% beta[-1L])
    beta[1] <- uniroot(function(a) sum(integration$w * plogis(a + slope)) - .375,
                       c(-20, 20), tol = 1e-12)$root
  }
  r <- plogis(as.vector(b %*% beta))
  continuous <- setting == "1trial_continuous_add"
  arm_means <- b %*% t(theta)
  if (!continuous) arm_means <- plogis(arm_means)
  truth <- sum(integration$w * (1 - r) * (arm_means[, 1] - arm_means[, 2])) /
    sum(integration$w * (1 - r))
  pr_s <- sim_softmax(b %*% t(gamma))
  expected_n <- colSums(pr_s * (row$N * integration$w * r))
  list(setting = setting, scenario = scenario, row = row, N = row$N,
       M = nrow(gamma), eta = eta, beta = beta, gamma = gamma, theta = theta,
       outcome = if (continuous) "continuous" else "binary", sigma = 1,
       truth = truth, expected_trial_n = expected_n, quadrature_order = order)
}

sim_softmax <- function(lp) {
  p <- exp(lp - apply(lp, 1L, max))
  p / rowSums(p)
}

sim_generate <- function(p) {
  u <- rnorm(p$N)
  z <- sqrt(p$eta[4]) * u + sqrt(1 - p$eta[4]) * matrix(rnorm(3L * p$N), p$N, 3L)
  d <- data.frame(X1 = as.integer(z[, 1] > qnorm(1 - p$eta[1])),
                  X2 = as.integer(z[, 2] > qnorm(1 - p$eta[2])), X3 = p$eta[3] + z[, 3])
  b <- sim_basis(d)
  d$R <- rbinom(p$N, 1L, plogis(as.vector(b %*% p$beta)))
  ii <- which(d$R == 1L)
  d$S <- 0L
  if (p$M == 1L) d$S[ii] <- 1L else {
    pr <- sim_softmax(b[ii, , drop = FALSE] %*% t(p$gamma))
    cumulative <- pr
    for (k in 2:p$M) cumulative[, k] <- cumulative[, k - 1L] + pr[, k]
    d$S[ii] <- 1L + rowSums(runif(length(ii)) > cumulative[, -p$M, drop = FALSE])
  }
  # Generate only observed trial outcomes; target outcomes are not required.
  d$A <- d$Y <- rep(NA_real_, p$N)
  d$A[ii] <- rbinom(length(ii), 1L, .5)
  means <- b[ii, , drop = FALSE] %*% t(p$theta)
  m <- means[cbind(seq_along(ii), ifelse(d$A[ii] == 1L, 1L, 2L))]
  d$Y[ii] <- if (p$outcome == "binary") rbinom(length(ii), 1L, plogis(m)) else
    m + rnorm(length(ii), sd = p$sigma)
  d
}

sim_aggregate <- function(d, p) {
  reports <- mu <- mu_se <- vector("list", p$M)
  names(mu) <- names(mu_se) <- paste0("S", seq_len(p$M))
  sizes <- integer(p$M)
  for (s in seq_len(p$M)) {
    x <- d[d$S == s, , drop = FALSE]
    sizes[s] <- nrow(x)
    if (nrow(x) < 2L) stop("Empty or singleton trial S", s)
    ind <- cbind(TRUE, x$X1 == 0, x$X1 == 1, x$X2 == 0, x$X2 == 1,
                 x$X3 <= 0, x$X3 > 0)
    ag <- data.frame(trial = paste0("S", s),
      subgroup = c("overall", "X1", "X1", "X2", "X2", "X3", "X3"),
      stratum = c("TRUE", "x == 0", "x == 1", "x == 0", "x == 1", "x <= 0", "x > 0"))
    for (a in 1:2) {
      ia <- ind & (x$A == if (a == 1L) 1L else 0L)
      na <- colSums(ia)
      ag[[paste0("n", a)]] <- na
      if (p$outcome == "binary") ag[[paste0("event", a)]] <- colSums(ia * x$Y) else {
        sy <- colSums(ia * x$Y)
        ag[[paste0("mean", a)]] <- sy / na
        ag[[paste0("se", a)]] <- sqrt(pmax(0, colSums(ia * x$Y^2) - sy^2 / na) /
                                        ((na - 1) * na))
      }
    }
    minimum <- if (p$outcome == "binary") 1L else 2L
    # Omit unreportable subgroups, but never silently omit an overall trial.
    keep <- ag$n1 >= minimum & ag$n2 >= minimum
    if (!keep[1L]) stop("Insufficient observations in an overall trial arm")
    reports[[s]] <- ag[keep, ]
    xbar <- colMeans(x[, c("X1", "X2", "X3")])
    sd3 <- sd(x$X3)
    mu[[s]] <- c(xbar, `I(X3^2)` = unname(xbar[3]^2 + (nrow(x) - 1) / nrow(x) * sd3^2))
    mu_se[[s]] <- c(X1 = unname(sqrt(xbar[1] * (1 - xbar[1]) / nrow(x))),
                   X2 = unname(sqrt(xbar[2] * (1 - xbar[2]) / nrow(x))),
                   X3 = sd3 / sqrt(nrow(x)), `I(X3^2)` = NA_real_)
  }
  list(reports = do.call(rbind, reports), mu = mu, mu_se = mu_se,
       trial_info = data.frame(trial = names(mu), n = sizes, allocation = .5))
}

sim_interval <- function(estimate, se) {
  if (length(estimate) != 1L || !is.finite(estimate)) stop("Nonfinite estimate")
  if (length(se) != 1L || !is.finite(se) || se < 0) stop("Invalid standard error")
  c(estimate = unname(estimate), se = unname(se),
    lower = unname(estimate - qnorm(.975) * se), upper = unname(estimate + qnorm(.975) * se))
}

sim_cima <- function(ag, target, outcome, return_fit = FALSE) {
  fit <- CIMAgD::CIMA(data_agd = ag$reports, mu = ag$mu, mu_se = ag$mu_se,
    Q = target, target = target, formula_cate = ~ X1 + X2 + X3,
    formula_tilting = ~ X1 + X2 + X3 + I(X3^2), outcome = outcome,
    trial_info = ag$trial_info, Sigma = NULL,
    q_target_shared = TRUE, jacobian_diagnostics = FALSE, tilting_diagnostics = FALSE,
    singularity_method = "none", match_subgroup_proportions = FALSE,
    time_indexed = FALSE, target_subgroups = list(), summaryX = NULL,
    q_method = "plugin", gmm_method = "one-step", variance = TRUE)
  # variance=TRUE is essential for coverage; bootstrap/copula fitting is unused.
  if (return_fit) return(fit)
  a <- fit$target$estimates
  a <- a[a$subgroup == "overall", ]
  sim_interval(a$estimate, a$se)
}

sim_ipd <- function(trial, target, outcome, pseudo = FALSE) {
  b <- sim_basis(trial); bt <- sim_basis(target)
  n <- nrow(b); nt <- nrow(bt)
  if (nt < 2L) stop("Insufficient target sample")
  if (pseudo) {
    y <- 2 * (2 * trial$A - 1) * trial$Y
    fit <- lm.fit(b, y)
    if (fit$rank < ncol(b)) stop("Singular pseudo-outcome regression")
    bread <- chol2inv(chol(crossprod(b)))
    # HC1 sandwich for the heteroskedastic pseudo-outcome, computed once.
    scores <- b * fit$residuals
    v <- n / (n - ncol(b)) * bread %*% crossprod(scores) %*% bread
    g <- as.vector(bt %*% fit$coefficients)
    gradient <- colMeans(bt)
  } else {
    design <- cbind(b, b * trial$A)
    if (outcome == "binary") {
      fit <- glm.fit(design, trial$Y, family = binomial(), control = glm.control(maxit = 50L))
      if (!fit$converged || fit$boundary || fit$rank < ncol(design) || any(!is.finite(fit$coefficients)))
        stop("Logistic g-formula did not yield a regular converged fit")
      coef <- fit$coefficients
      p0 <- plogis(as.vector(bt %*% coef[1:4]))
      p1 <- plogis(as.vector(bt %*% (coef[1:4] + coef[5:8])))
      g <- p1 - p0
      gradient <- c(colMeans(bt * (p1 * (1 - p1) - p0 * (1 - p0))),
                    colMeans(bt * (p1 * (1 - p1))))
      v <- chol2inv(chol(crossprod(design, design * (fit$fitted.values * (1 - fit$fitted.values)))))
    } else {
      fit <- lm.fit(design, trial$Y)
      if (fit$rank < ncol(design)) stop("Singular continuous g-formula regression")
      v <- sum(fit$residuals^2) / (n - ncol(design)) * chol2inv(chol(crossprod(design)))
      g <- as.vector(bt %*% fit$coefficients[5:8])
      gradient <- c(rep(0, 4), colMeans(bt))
    }
  }
  sim_interval(mean(g), sqrt(as.numeric(crossprod(gradient, v %*% gradient)) + var(g) / nt))
}

sim_meta <- function(ag, target, regression = FALSE) {
  d <- ag$reports[ag$reports$subgroup == "overall", ]
  yi <- d$event1 / d$n1 - d$event2 / d$n2
  vi <- (d$event1 / d$n1) * (1 - d$event1 / d$n1) / d$n1 +
    (d$event2 / d$n2) * (1 - d$event2 / d$n2) / d$n2
  if (any(!is.finite(vi) | vi <= 0)) stop("Nonpositive overall effect variance")
  if (!regression) {
    fit <- meta::metagen(TE = yi, seTE = sqrt(vi), studlab = d$trial, sm = "RD",
      common = FALSE, random = TRUE, method.tau = "REML", method.tau.ci = "",
      method.random.ci = "classic", prediction = FALSE, overall.hetstat = FALSE,
      keepdata = FALSE, keeprma = FALSE, backtransf = FALSE)
    return(sim_interval(fit$TE.random, fit$seTE.random))
  }
  mods <- t(vapply(ag$mu[d$trial], function(x) x[c("X1", "X2", "X3")], numeric(3)))
  if (qr(cbind(1, mods))$rank < 4L) stop("Singular trial-level meta-regression")
  fit <- metafor::rma.uni(yi = yi, vi = vi, mods = mods, method = "REML", test = "z",
                         skipr2 = TRUE, outlist = "minimal")
  x0 <- c(1, colMeans(target))
  vtarget <- as.numeric(crossprod(fit$b[-1], cov(target) %*% fit$b[-1])) / nrow(target)
  sim_interval(sum(x0 * fit$b), sqrt(as.numeric(crossprod(x0, fit$vb %*% x0)) + vtarget))
}

sim_capture <- function(expr) {
  warnings <- character()
  error <- ""
  result <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    c(estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_)
  })
  data.frame(as.list(result), success = !nzchar(error), error = error,
             warning = paste(unique(warnings), collapse = " | "), stringsAsFactors = FALSE)
}

sim_replication <- function(replication, p, seed) {
  # Disjoint deterministic seeds: same draw for a replication regardless of
  # worker count, selected scenario order, or splitting a run into batches.
  rng_seed <- as.integer((as.double(seed) + match(p$setting, sim_settings) * 10000000 +
                          p$scenario * 100000 + replication) %% 2147483646 + 1)
  set.seed(rng_seed, kind = "L'Ecuyer-CMRG")
  methods <- c("CIMAgD", "IPD-PO", "IPD-GF", if (p$M == 5L) c("Meta", "Meta-reg"))
  t0 <- proc.time()[[3]]
  d <- sim_generate(p)
  trial <- d[d$R == 1L, ]; target <- d[d$R == 0L, c("X1", "X2", "X3")]
  ag <- tryCatch(sim_aggregate(d, p), error = function(e) e)
  need_ag <- function() { if (inherits(ag, "error")) stop(conditionMessage(ag)); ag }
  result <- list(sim_capture(sim_cima(need_ag(), target, p$outcome)),
                 sim_capture(sim_ipd(trial, target, p$outcome, pseudo = TRUE)),
                 sim_capture(sim_ipd(trial, target, p$outcome)))
  if (p$M == 5L) result <- c(result, list(sim_capture(sim_meta(need_ag(), target)),
                                         sim_capture(sim_meta(need_ag(), target, regression = TRUE))))
  result <- do.call(rbind, result)
  sizes <- tabulate(d$S[d$R == 1L], nbins = p$M)
  cbind(setting = p$setting, scenario = p$scenario, replication = replication,
        seed = rng_seed, method = methods, truth = p$truth, result,
        n_target = nrow(target), n_trial = nrow(trial),
        min_trial_n = min(sizes), max_trial_n = max(sizes),
        n_reports = if (inherits(ag, "error")) NA_integer_ else nrow(ag$reports),
        elapsed_seconds = proc.time()[[3]] - t0, row.names = NULL)
}

sim_check_packages <- function(setting) {
  packages <- c("CIMAgD", if (startsWith(setting, "5trial")) c("meta", "metafor"))
  for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE))
    stop("Install package '", pkg, "' before running this setting.")
  required <- c("data_agd", "formula_cate", "formula_tilting", "mu_se", "q_target_shared",
                "jacobian_diagnostics", "tilting_diagnostics", "gmm_method", "variance")
  if (!all(required %in% names(formals(CIMAgD::CIMA))))
    stop("Installed CIMAgD has an old interface; install the updated local package first.")
  setNames(vapply(packages, function(x) as.character(utils::packageVersion(x)), ""), packages)
}

sim_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("--help", "-h"))) {
    cat("Default: 2 replications of scenarios 1 and 16 (smoke only).\n",
        "--full                 1000 replications of all 16 scenarios\n",
        "--reps=INTEGER         Replications per selected scenario (1..100000)\n",
        "--scenarios=1,2,16     Explicit Appendix scenario indices (or all)\n",
        "--workers=INTEGER      PSOCK workers (default 1; Windows/macOS/Linux)\n",
        "--seed=INTEGER         Master seed (default 20261019)\n",
        "--out=PATH             Root output directory\n",
        "--overwrite            Replace matching scenario output files\n", sep = "")
    return(NULL)
  }
  full <- "--full" %in% args
  opts <- list(reps = if (full) 1000L else 2L, scenarios = if (full) 1:16 else c(1L, 16L),
               workers = 1L, seed = 20261019L, out = if (full) "results" else "results_smoke",
               overwrite = "--overwrite" %in% args)
  for (a in setdiff(args, c("--full", "--overwrite"))) {
    if (!grepl("^--[^=]+=.+$", a)) stop("Unrecognized argument: ", a)
    key <- sub("^--([^=]+)=.*$", "\\1", a); value <- sub("^--[^=]+=", "", a)
    if (!key %in% c("reps", "scenarios", "workers", "seed", "out")) stop("Unknown option: ", key)
    opts[[key]] <- switch(key, out = value,
      scenarios = if (value == "all") 1:16 else as.numeric(strsplit(value, ",", fixed = TRUE)[[1]]),
      as.numeric(value))
  }
  for (key in c("reps", "workers", "seed")) {
    x <- opts[[key]]
    if (length(x) != 1L || !is.finite(x) || x != floor(x) || x < 1 || x > .Machine$integer.max)
      stop("Invalid --", key)
  }
  if (opts$reps > 100000L) stop("--reps must not exceed 100000")
  if (!length(opts$scenarios) || anyNA(opts$scenarios) || any(!opts$scenarios %in% 1:16) ||
      anyDuplicated(opts$scenarios)) stop("Scenario indices must be unique integers in 1..16")
  opts
}

sim_run <- function(setting, script_dir, args = commandArgs(trailingOnly = TRUE)) {
  opts <- sim_options(args)
  if (is.null(opts)) return(invisible(NULL))
  versions <- sim_check_packages(setting)
  implementation_md5 <- unname(tools::md5sum(file.path(script_dir, "simulation_common.R")))
  out <- file.path(opts$out, setting)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(out, sprintf("scenario_%02d.rds", opts$scenarios))
  if (any(file.exists(files)) && !opts$overwrite)
    stop("Output already exists; use another --out directory or explicitly --overwrite.")
  workers <- min(opts$workers, opts$reps)
  cl <- NULL
  if (workers > 1L) {
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(path, libraries, setting) {
      .libPaths(libraries); sys.source(path, envir = .GlobalEnv)
      sim_check_packages(setting); NULL
    }, normalizePath(file.path(script_dir, "simulation_common.R")), .libPaths(), setting)
  }
  for (index in seq_along(opts$scenarios)) {
    s <- opts$scenarios[index]
    p <- sim_parameters(setting, s)
    t0 <- proc.time()[[3]]
    run <- if (is.null(cl)) lapply(seq_len(opts$reps), sim_replication, p = p, seed = opts$seed) else
      parallel::parLapply(cl, seq_len(opts$reps), sim_replication, p = p, seed = opts$seed)
    results <- do.call(rbind, run)
    payload <- list(schema_version = 1L, setting = setting, scenario_table = sim_grid(setting),
                    parameters = p, options = opts, package_versions = versions,
                    implementation_md5 = implementation_md5,
                    session_info = utils::sessionInfo(), results = results)
    # One checkpoint per scenario; never retain the large fitted objects.
    temp <- tempfile(pattern = "scenario_", tmpdir = out)
    saveRDS(payload, temp, compress = FALSE)
    if (!file.rename(temp, files[index])) {
      if (!file.copy(temp, files[index], overwrite = opts$overwrite)) stop("Cannot save ", files[index])
      unlink(temp)
    }
    failures <- sum(!results$success)
    message(sprintf("%s: scenario %02d, %d replications, %d/%d failed method fits, %.1f seconds",
                    setting, s, opts$reps, failures, nrow(results), proc.time()[[3]] - t0))
  }
  invisible(files)
}
