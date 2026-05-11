# validate_irt_2pl.R
# Run from library root after devtools::load_all()
# or: source("R/irt_2pl.R"); source("R/iDIFr-package.R")

library(TAM)

cat("=============================================================\n")
cat("Validation 1: Single-group parameter recovery\n")
cat("=============================================================\n\n")

set.seed(42)
dat <- simulate_dif(n_persons=1000, n_items=10, n_groups=1,
                    dif_items=c(), dif_effect=0, seed=42)
resp <- as.matrix(dat[, paste0("item_", 1:10)])

cat("Fitting iDIFr 2PL...\n")
t0 <- proc.time()
m_i <- fit_2pl(resp, verbose=FALSE)
cat(sprintf("  Time: %.2fs  |  iter: %d  |  converged: %s  |  loglik: %.4f\n\n",
            (proc.time()-t0)["elapsed"], m_i$iter, m_i$converged, m_i$loglik))

cat("Fitting TAM 2PL...\n")
t0 <- proc.time()
m_t <- suppressMessages(TAM::tam.mml.2pl(resp=resp, verbose=FALSE))
cat(sprintf("  Time: %.2fs  |  loglik: %.4f\n\n",
            (proc.time()-t0)["elapsed"], m_t$ic$loglike))

cat("Comparison:\n")
cmp <- data.frame(
  item   = paste0("item_", 1:10),
  a_i    = round(m_i$item_params$a, 4),
  a_t    = round(m_t$item_irt$alpha, 4),
  a_diff = round(m_i$item_params$a - m_t$item_irt$alpha, 4),
  b_i    = round(m_i$item_params$b, 4),
  b_t    = round(m_t$item_irt$beta, 4),
  b_diff = round(m_i$item_params$b - m_t$item_irt$beta, 4)
)
print(cmp)
cat(sprintf("\nMAD a: %.4f  |  MAD b: %.4f\n",
            mean(abs(cmp$a_diff)), mean(abs(cmp$b_diff))))
cat(sprintf("loglik diff: %.4f  (iDIFr - TAM)\n\n",
            m_i$loglik - m_t$ic$loglike))

cat("=============================================================\n")
cat("Validation 2: Multigroup constrained recovery\n")
cat("=============================================================\n\n")

set.seed(42)
dat2  <- simulate_dif(n_persons=1000, n_items=10, n_groups=2,
                      dif_items=c(), dif_effect=0, seed=42)
resp2  <- as.matrix(dat2[, paste0("item_", 1:10)])
grp2   <- dat2$group

cat("Fitting iDIFr 2PL (constrained)...\n")
t0 <- proc.time()
m_i2 <- fit_2pl(resp2, group=grp2, constrain="items", verbose=FALSE)
cat(sprintf("  Time: %.2fs  |  iter: %d  |  converged: %s  |  loglik: %.4f\n\n",
            (proc.time()-t0)["elapsed"], m_i2$iter, m_i2$converged, m_i2$loglik))

cat("Fitting TAM 2PL (multigroup)...\n")
t0 <- proc.time()
m_t2 <- suppressMessages(TAM::tam.mml.2pl(resp=resp2,
                                           group=factor(grp2), verbose=FALSE))
cat(sprintf("  Time: %.2fs  |  loglik: %.4f\n\n",
            (proc.time()-t0)["elapsed"], m_t2$ic$loglike))

cat("Comparison:\n")
cmp2 <- data.frame(
  item   = paste0("item_", 1:10),
  a_i    = round(m_i2$item_params$a, 4),
  a_t    = round(m_t2$item_irt$alpha, 4),
  a_diff = round(m_i2$item_params$a - m_t2$item_irt$alpha, 4),
  b_i    = round(m_i2$item_params$b, 4),
  b_t    = round(m_t2$item_irt$beta, 4),
  b_diff = round(m_i2$item_params$b - m_t2$item_irt$beta, 4)
)
print(cmp2)
cat(sprintf("\nMAD a: %.4f  |  MAD b: %.4f\n",
            mean(abs(cmp2$a_diff)), mean(abs(cmp2$b_diff))))
cat(sprintf("loglik diff: %.4f  (iDIFr - TAM)\n", m_i2$loglik - m_t2$ic$loglike))
cat(sprintf("Group params: G1 mu=%.4f s2=%.4f  |  G2 mu=%.4f s2=%.4f\n\n",
            m_i2$group_params$mu[1], m_i2$group_params$sigma2[1],
            m_i2$group_params$mu[2], m_i2$group_params$sigma2[2]))

cat("=============================================================\n")
cat("Validation 3: DIF detection — items 3 and 7 have DIF\n")
cat("=============================================================\n\n")

set.seed(42)
dat3 <- simulate_dif(n_persons=1000, n_items=10, n_groups=2,
                     dif_items=c(3,7), dif_effect=1.2, seed=42)
resp3 <- as.matrix(dat3[, paste0("item_", 1:10)])
grp3  <- dat3$group

cat("Fitting constrained model...\n")
m_c <- fit_2pl(resp3, group=grp3, constrain="items",  verbose=FALSE)
cat(sprintf("  loglik: %.4f\n", m_c$loglik))

cat("Fitting free model (constrain='none')...\n")
m_f <- fit_2pl(resp3, group=grp3, constrain="none", verbose=FALSE)
cat(sprintf("  loglik: %.4f\n\n", m_f$loglik))

# Per-item LL
ll_c <- item_loglik_mg(m_c)
ll_f <- numeric(10)
for (gi in 1:2) {
  idx    <- which(grp3 == m_f$group_levels[gi])
  resp_g <- resp3[idx, ]
  post_g <- m_f$posterior[idx, ]
  ll_f   <- ll_f + item_loglik(m_f, resp_g, post_g, gi)
}

chi_sq <- pmax(0, -2*(ll_c - ll_f))
df     <- 2*(2-1)
p_val  <- pchisq(chi_sq, df, lower.tail=FALSE)

cat("Per-item LRT (items 3 and 7 should be significant):\n")
res <- data.frame(
  item   = paste0("item_", 1:10),
  chi_sq = round(chi_sq, 3),
  df     = df,
  p_value = round(p_val, 4),
  flag   = ifelse(p_val < 0.05, "DIF *", "ok")
)
print(res)

cat("\n=============================================================\n")
cat("Validation 4: Uniform vs non-uniform DIF\n")
cat("=============================================================\n\n")

cat("Fitting alpha-constrained model (uniform DIF test)...\n")
m_alpha <- fit_2pl(resp3, group=grp3, constrain="alpha", verbose=FALSE)
cat(sprintf("  loglik: %.4f\n", m_alpha$loglik))

cat("Fitting beta-constrained model (non-uniform DIF test)...\n")
m_beta <- fit_2pl(resp3, group=grp3, constrain="beta", verbose=FALSE)
cat(sprintf("  loglik: %.4f\n\n", m_beta$loglik))

# Uniform DIF:     constrained vs alpha-constrained (b freed)
# Non-uniform DIF: alpha-constrained vs fully free (a freed beyond uniform)
ll_alpha_mg <- item_loglik_mg(m_alpha)
ll_beta_mg  <- item_loglik_mg(m_beta)

ll_alpha_f <- numeric(10)
ll_beta_f  <- numeric(10)
for (gi in 1:2) {
  idx    <- which(grp3 == m_f$group_levels[gi])
  resp_g <- resp3[idx, ]
  post_g <- m_f$posterior[idx, ]
  ll_alpha_f <- ll_alpha_f + item_loglik(m_f, resp_g, post_g, gi)
  ll_beta_f  <- ll_beta_f  + item_loglik(m_f, resp_g, post_g, gi)
}

chi_uniform    <- pmax(0, -2*(ll_c - ll_alpha_f))
chi_nonuniform <- pmax(0, -2*(ll_alpha_f - ll_f))
p_uniform      <- pchisq(chi_uniform,    df=1, lower.tail=FALSE)
p_nonuniform   <- pchisq(chi_nonuniform, df=1, lower.tail=FALSE)

cat("Uniform vs non-uniform decomposition:\n")
res2 <- data.frame(
  item      = paste0("item_", 1:10),
  chi_unif  = round(chi_uniform, 3),
  p_unif    = round(p_uniform, 4),
  chi_nonunif = round(chi_nonuniform, 3),
  p_nonunif = round(p_nonuniform, 4)
)
print(res2)
cat("\n(Items generated with uniform DIF — chi_unif should dominate)\n")
