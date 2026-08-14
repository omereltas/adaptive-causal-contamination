# =============================================================================
# ADAPTIVE MONTE CARLO FRAMEWORK - PRODUCTION RUN (n_mc = 2000)
# =============================================================================
# BU SURUME KADAR UYGULANAN TUM DUZELTMELER:
#   [FIX 1] beta_ps / beta_out vektorleri 11 -> 21 elemana genisletildi.
#           (Eskiden p=20 icin sessizce NA donuyordu, tum p=20 senaryolari
#            hic yazilmadan atlaniyordu - bu duzeltilmeden ONCE kesfedildi
#            ve ayri bir kosuyla (run_p20_fix.R) telafi edildi.)
#   [FIX 2] select_estimator_adaptive(): covMcd() -> guvenli IQR-tabanli
#           outlier tespiti. moments::kurtosis() eksikse manuel fallback.
#   [FIX 3] Checkpoint dosyasi olusturma bug'i giderildi: eskiden
#           write.csv(data.frame(), file) ile olusturulan "bos" stub
#           dosyasi aslinda 1 satir icerdigi (file.size>0) icin ilk
#           gercek yazimda header HIC YAZILMIYORDU. Artik sadece
#           file.exists() kontrolu kullaniliyor, stub dosyasi hic
#           olusturulmuyor.
#   [FIX 4] n_mc TUTARLILIK KONTROLU eklendi: eger checkpoint dosyasi
#           zaten varsa ve icindeki n_mc degeri CONFIG$n_mc ile
#           UYUSMUYORSA script DURDURULUR (eskiden sessizce n_mc=10'luk
#           "tamamlanmis" senaryolari atlayip yanlislikla dusuk kaliteli
#           sonuclari "nihai" sonuc olarak birakma riski vardi).
#   [FIX 5] doRNG eklendi: paralel worker'larin RNG akislari artik
#           bagimsiz ve tekrarlanabilir (L'Ecuyer-CMRG). Eskiden %dopar%
#           ile worker'lar arasi RNG bagimsizligi GARANTI DEGILDI - JASA/
#           Biometrika hakemleri bunu sorabilir, simdi guvenli.
#
# BILINCLI OLARAK YAPILMAYAN DEGISIKLIK:
#   Senaryo x tekrar duzlestirme (onceden konusulan "1 numarali duzeltme")
#   UYGULANMADI. Gerekce: n_mc=2000/3 cekirdek=666 gorev/cekirdek zaten
#   iyi dengeleniyor (n_mc=10 iken sorun olan yuk dengesizligi artik
#   onemsiz). 675 senaryonun ayri ayri dispatch edilmesinin ek maliyeti
#   (birkac saniye) ~10 gunluk toplam surenin yaninda ihmal edilebilir.
#   Test edilmis/kanitlanmis mevcut checkpoint granularitesini (senaryo
#   bazinda) bozup riske atmaya gerek yok.
#
# TAHMINI SURE: ~9.9 gun (675 senaryo, gercek n_mc=10 olcumunden 200x
# olceklenerek, 3 cekirdek - Intel i5-4460). KAPSAM DARALTILMADI.
#
# CALISTIRMA ONERISI: Bu scripti RStudio icinde DEGIL, komut satirindan
# calistirin (asagida "CALISTIRMA TALIMATLARI" bolumune bakin) - gunlerce
# surecek bir isi interaktif oturuma bagli tutmak riskli.
# =============================================================================

rm(list = ls())

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
CONFIG <- list(
  n_mc = 2000,                        # HEDEF: tam replikasyon sayisi
  n_values = c(200, 500, 1000),
  p_values = c(5, 10, 20),
  true_ate = 2.0,
  contamination_types = c("none", "covariate", "outcome", "treatment", "mixed"),
  contamination_rates = c(0, 0.05, 0.10, 0.15, 0.20),
  contamination_strengths = c(3, 5, 10),
  n_cores = 3,                        # i5-4460: 4 fiziksel cekirdek, HT yok -> detectCores()-1=3
  # ONEMLI: n_mc=10 test kosularindan FARKLI dosya adlari - eski checkpoint
  # ile KARISMASIN diye. Bu isimler zaten "production" icin ayrilmis olsun.
  checkpoint_file = "mc_checkpoint_production_nmc2000.csv",
  results_file    = "mc_results_production_nmc2000.csv",
  seed = 42,
  verbose = TRUE
)

# =============================================================================
# 1. PACKAGES
# =============================================================================
packages <- c("MASS", "robustbase", "sandwich", "lmtest", "ggplot2",
              "reshape2", "gridExtra", "parallel", "doParallel", "foreach",
              "doRNG", "moments", "MatchIt")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

USE_DORNG <- requireNamespace("doRNG", quietly = TRUE)
if (!USE_DORNG) {
  warning("doRNG yuklenemedi - paralel RNG akislarinin bagimsizligi GARANTI DEGIL. ",
          "Sonuclarin tekrarlanabilirligi etkilenebilir.")
}

set.seed(CONFIG$seed)

# =============================================================================
# 2. THEORETICAL SECTION: IF / ASYMPTOTIC BIAS / BREAKDOWN POINT
#    (degismedi - orijinal kod, sadece generate_data/calc_* asagida FIXED)
# =============================================================================
calc_if_ipw <- function(data, perturbation_amount = NULL) {
  n <- nrow(data)
  if (is.null(perturbation_amount)) perturbation_amount <- 1/n
  ate_full <- calc_ipw_naive(data)
  if_values <- numeric(n)
  for (i in 1:n) {
    dp <- data
    dp$Y[i] <- dp$Y[i] + perturbation_amount * sd(data$Y)
    tryCatch({
      if_values[i] <- (calc_ipw_naive(dp) - ate_full) / perturbation_amount
    }, error = function(e) { if_values[i] <<- 0 })
  }
  if_values
}

calc_if_ipw_robust <- function(data, perturbation_amount = NULL) {
  n <- nrow(data)
  if (is.null(perturbation_amount)) perturbation_amount <- 1/n
  ate_full <- calc_ipw_robust(data)
  if_values <- numeric(n)
  for (i in 1:n) {
    dp <- data
    dp$Y[i] <- dp$Y[i] + perturbation_amount * sd(data$Y)
    tryCatch({
      if_values[i] <- (calc_ipw_robust(dp) - ate_full) / perturbation_amount
    }, error = function(e) { if_values[i] <<- 0 })
  }
  if_values
}

calc_if_dr <- function(data, perturbation_amount = NULL) {
  n <- nrow(data)
  if (is.null(perturbation_amount)) perturbation_amount <- 1/n
  ate_full <- calc_dr_robust(data)
  if_values <- numeric(n)
  for (i in 1:n) {
    dp <- data
    dp$Y[i] <- dp$Y[i] + perturbation_amount * sd(data$Y)
    tryCatch({
      if_values[i] <- (calc_dr_robust(dp) - ate_full) / perturbation_amount
    }, error = function(e) { if_values[i] <<- 0 })
  }
  if_values
}

plot_influence_functions <- function(data) {
  cat("Computing Influence Functions...\n")
  if_ipw <- calc_if_ipw(data); if_ipw_r <- calc_if_ipw_robust(data); if_dr <- calc_if_dr(data)
  if_df <- data.frame(Observation = 1:nrow(data), IPW = if_ipw, IPW_Robust = if_ipw_r, DR_Robust = if_dr)
  if_df_long <- reshape2::melt(if_df, id.vars = "Observation")
  p <- ggplot(if_df_long, aes(x = Observation, y = value, color = variable)) +
    geom_line(alpha = 0.7) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Empirical Influence Functions", x = "Observation Index", y = "Influence Value", color = "Estimator") +
    theme_minimal()
  cat(sprintf("Max |IF| naive=%.3f robust_ipw=%.3f dr=%.3f\n", max(abs(if_ipw)), max(abs(if_ipw_r)), max(abs(if_dr))))
  list(plot = p, if_data = if_df)
}

calc_asymptotic_bias <- function(data_clean, data_contaminated, epsilon) {
  true_ate <- calc_ipw_naive(data_clean)
  bias_naive <- calc_ipw_naive(data_contaminated) - true_ate
  bias_ipw_r <- calc_ipw_robust(data_contaminated) - true_ate
  bias_dr_r  <- calc_dr_robust(data_contaminated) - true_ate
  bias_aipw_r<- calc_aipw_robust(data_contaminated) - true_ate
  if_naive <- calc_if_ipw(data_contaminated); if_ipw_r <- calc_if_ipw_robust(data_contaminated)
  n <- nrow(data_contaminated); n_contam <- round(epsilon * n)
  contam_idx <- order(abs(if_naive), decreasing = TRUE)[1:n_contam]
  result <- data.frame(
    Method = c("Naive IPW", "Robust IPW", "DR Robust", "AIPW Robust"),
    Observed_Bias = c(bias_naive, bias_ipw_r, bias_dr_r, bias_aipw_r),
    Theoretic_Bias_1stOrder = c(epsilon * mean(if_naive[contam_idx]), epsilon * mean(if_ipw_r[contam_idx]), NA, NA),
    Epsilon = epsilon
  )
  print(result)
  result
}

calc_breakdown_point <- function(data, method = "ipw_naive", max_epsilon = 0.5, step = 0.02, bias_threshold = NULL) {
  true_ate <- calc_ipw_naive(data)
  if (is.null(bias_threshold)) bias_threshold <- 2 * abs(true_ate)
  epsilons <- seq(0, max_epsilon, by = step)
  biases <- numeric(length(epsilons))
  for (i in seq_along(epsilons)) {
    sim_data <- generate_data(n = nrow(data), p = ncol(data) - 2, true_ate = true_ate,
                              contamination_type = "mixed", contamination_rate = epsilons[i], contamination_strength = 10)
    est <- tryCatch(switch(method, "ipw_naive" = calc_ipw_naive(sim_data$data), "ipw_robust" = calc_ipw_robust(sim_data$data),
                           "dr_robust" = calc_dr_robust(sim_data$data), "aipw_robust" = calc_aipw_robust(sim_data$data), NA),
                    error = function(e) NA)
    biases[i] <- ifelse(is.na(est), Inf, abs(est - true_ate))
    if (biases[i] > bias_threshold) break
  }
  bp <- epsilons[max(which(biases[1:i] <= bias_threshold))]   # [FIX] sadece doldurulmus kismi kullan
  cat(sprintf("BP(%s) = %.2f\n", method, bp))
  list(breakdown_point = bp, epsilons = epsilons[1:i], biases = biases[1:i])
}

compare_breakdown_points <- function(data) {
  methods <- c("ipw_naive", "ipw_robust", "dr_robust", "aipw_robust")
  bp_results <- lapply(methods, function(m) calc_breakdown_point(data, method = m))
  names(bp_results) <- methods
  bp_df <- do.call(rbind, lapply(names(bp_results), function(m)
    data.frame(Method = m, Epsilon = bp_results[[m]]$epsilons, Bias = bp_results[[m]]$biases)))
  p <- ggplot(bp_df, aes(x = Epsilon, y = Bias, color = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2) +
    geom_hline(yintercept = 2 * abs(calc_ipw_naive(data)), linetype = "dashed", color = "red") +
    labs(title = "Breakdown Point Analysis", x = "Contamination Rate", y = "|Bias|") + theme_minimal()
  list(plot = p, bp_results = bp_results)
}

# =============================================================================
# 3. DATA GENERATION - [FIX 1] BETA VEKTORLERI GENISLETILDI (p=20 guvenli)
# =============================================================================
beta_ps_base  <- c(0.5, -0.3, 0.4, -0.2, 0.1, 0.3, -0.1, 0.2, -0.15, 0.1, 0.05)
beta_ps_extra <- c(-0.08, 0.06, -0.05, 0.04, -0.03, 0.03, -0.02, 0.02, -0.01, 0.01)
BETA_PS_FULL  <- c(beta_ps_base, beta_ps_extra)   # uzunluk 21

beta_out_base  <- c(1.0, 0.5, -0.4, 0.3, -0.2, 0.1, -0.05, 0.08, -0.03, 0.06, 0.02)
beta_out_extra <- c(0.04, -0.03, 0.025, -0.02, 0.015, -0.01, 0.01, -0.008, 0.006, -0.005)
BETA_OUT_FULL  <- c(beta_out_base, beta_out_extra) # uzunluk 21

generate_data <- function(n, p = 5, true_ate = 2.0, contamination_type = "none",
                          contamination_rate = 0.1, contamination_strength = 5.0,
                          covariate_type = "normal") {
  stopifnot(p + 1 <= length(BETA_PS_FULL), p + 1 <= length(BETA_OUT_FULL))  # [FIX] sessiz NA yerine hata
  
  Sigma <- matrix(0.3, p, p); diag(Sigma) <- 1.0
  if (covariate_type == "normal") {
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  } else if (covariate_type == "t") {
    X <- rmvt(n, sigma = Sigma, df = 5)
  } else if (covariate_type == "skewed") {
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma); X <- exp(X) - 1
  }
  colnames(X) <- paste0("X", 1:p)
  
  beta_ps <- BETA_PS_FULL[1:(p+1)]
  logit_ps <- cbind(1, X) %*% beta_ps
  if (p >= 2) logit_ps <- logit_ps + 0.1 * X[,1] * X[,2]
  T <- rbinom(n, 1, plogis(logit_ps))
  
  beta_out <- BETA_OUT_FULL[1:(p+1)]
  Y0 <- cbind(1, X) %*% beta_out + rnorm(n, 0, 1)
  if (p >= 2) Y0 <- Y0 + 0.2 * X[,1]^2 - 0.1 * X[,2]^2
  Y1 <- Y0 + true_ate + rnorm(n, 0, 0.5)
  Y <- T * Y1 + (1 - T) * Y0
  
  n_contam <- floor(n * contamination_rate)
  if (n_contam > 0) {
    contam_idx <- sample(1:n, n_contam)
    if (contamination_type == "covariate") {
      X[contam_idx, ] <- X[contam_idx, ] + contamination_strength * matrix(rt(n_contam * p, df = 3), n_contam, p)
    } else if (contamination_type == "outcome") {
      Y[contam_idx] <- Y[contam_idx] + contamination_strength * rt(n_contam, df = 3)
    } else if (contamination_type == "treatment") {
      T[contam_idx] <- 1 - T[contam_idx]
    } else if (contamination_type == "mixed") {
      n_each <- max(1, floor(n_contam / 3))
      idx1 <- contam_idx[1:min(n_each, length(contam_idx))]
      idx2 <- contam_idx[(n_each+1):min(2*n_each, length(contam_idx))]
      idx3 <- contam_idx[(2*n_each+1):length(contam_idx)]
      if (length(idx1) > 0) X[idx1, ] <- X[idx1, ] + contamination_strength * matrix(rt(length(idx1) * p, df = 3), length(idx1), p)
      if (length(idx2) > 0) Y[idx2] <- Y[idx2] + contamination_strength * rt(length(idx2), df = 3)
      if (length(idx3) > 0) T[idx3] <- 1 - T[idx3]
    } else if (contamination_type == "leverage") {
      X[contam_idx, ] <- X[contam_idx, ] + contamination_strength * 2
      Y[contam_idx] <- Y[contam_idx] + contamination_strength * rnorm(n_contam, 0, 1)
    }
  }
  data <- data.frame(Y = as.vector(Y), T = T, X)
  list(data = data, true_ate = true_ate, contamination_idx = if (n_contam > 0) contam_idx else integer(0))
}

# =============================================================================
# 4. ESTIMATORS (degismedi)
# =============================================================================
calc_ipw_robust <- function(data, trim_lower = 0.01, trim_upper = 0.99, huber_k = 1.345) {
  ps_formula <- as.formula(paste("T ~", paste(names(data)[-(1:2)], collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps_hat <- pmax(trim_lower, pmin(trim_upper, predict(ps_model, type = "response")))
  w1 <- data$T / ps_hat; w0 <- (1 - data$T) / (1 - ps_hat)
  y1 <- data$Y[data$T == 1]; y0 <- data$Y[data$T == 0]
  med1 <- median(y1, na.rm = TRUE); med0 <- median(y0, na.rm = TRUE)
  mad1 <- max(mad(y1, na.rm = TRUE), 1e-10); mad0 <- max(mad(y0, na.rm = TRUE), 1e-10)
  hw <- function(r, k) pmin(1, k / abs(r))
  w1_r <- w1 * hw((data$Y - med1) / mad1, huber_k); w0_r <- w0 * hw((data$Y - med0) / mad0, huber_k)
  sum(w1_r * data$Y) / sum(w1_r) - sum(w0_r * data$Y) / sum(w0_r)
}

calc_dr_robust <- function(data, trim_lower = 0.01, trim_upper = 0.99) {
  ps_formula <- as.formula(paste("T ~", paste(names(data)[-(1:2)], collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps_hat <- pmax(trim_lower, pmin(trim_upper, predict(ps_model, type = "response")))
  cov_names <- names(data)[-(1:2)]
  out_model <- lmrob(as.formula(paste("Y ~", paste(c("T", cov_names), collapse = " + "))), data = data, setting = "KS2014")
  d1 <- data; d1$T <- 1; d0 <- data; d0$T <- 0
  mu1 <- predict(out_model, newdata = d1); mu0 <- predict(out_model, newdata = d0)
  r1 <- data$Y - mu1; r0 <- data$Y - mu0
  m1 <- max(mad(r1, na.rm = TRUE), 1e-10); m0 <- max(mad(r0, na.rm = TRUE), 1e-10)
  hw <- function(r, k = 1.345) pmin(1, k / abs(r))
  mean(mu1 + hw(r1/m1) * data$T * r1 / ps_hat) - mean(mu0 + hw(r0/m0) * (1 - data$T) * r0 / (1 - ps_hat))
}

calc_aipw_robust <- function(data, trim_lower = 0.01, trim_upper = 0.99) {
  ps_formula <- as.formula(paste("T ~", paste(names(data)[-(1:2)], collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps_hat <- pmax(trim_lower, pmin(trim_upper, predict(ps_model, type = "response")))
  cov_names <- names(data)[-(1:2)]
  of <- as.formula(paste("Y ~", paste(cov_names, collapse = " + ")))
  m1 <- lmrob(of, data = data[data$T == 1, ], setting = "KS2014")
  m0 <- lmrob(of, data = data[data$T == 0, ], setting = "KS2014")
  mu1 <- predict(m1, newdata = data); mu0 <- predict(m0, newdata = data)
  eif <- (data$T * (data$Y - mu1) / ps_hat + mu1) - ((1 - data$T) * (data$Y - mu0) / (1 - ps_hat) + mu0)
  median(eif)
}

calc_ipw_naive <- function(data, trim_lower = 0.01, trim_upper = 0.99) {
  ps_formula <- as.formula(paste("T ~", paste(names(data)[-(1:2)], collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps_hat <- pmax(trim_lower, pmin(trim_upper, predict(ps_model, type = "response")))
  w1 <- data$T / ps_hat; w0 <- (1 - data$T) / (1 - ps_hat)
  sum(w1 * data$Y) / sum(w1) - sum(w0 * data$Y) / sum(w0)
}

calc_reg_adj <- function(data) {
  cov_names <- names(data)[-(1:2)]
  model <- lm(as.formula(paste("Y ~", paste(c("T", cov_names), collapse = " + "))), data = data)
  d1 <- data; d1$T <- 1; d0 <- data; d0$T <- 0
  mean(predict(model, newdata = d1)) - mean(predict(model, newdata = d0))
}

# =============================================================================
# 5. ADAPTIVE SELECTION - [FIX 2] covMcd -> IQR, moments fallback
# =============================================================================
select_estimator_adaptive <- function(data) {
  X <- data[, -(1:2)]; p <- ncol(X)
  
  cov_outlier_rate <- 0
  tryCatch({
    outlier_flags <- rep(0, nrow(X))
    for (j in 1:ncol(X)) {
      xj <- X[, j]; q1 <- quantile(xj, 0.25, na.rm = TRUE); q3 <- quantile(xj, 0.75, na.rm = TRUE); iqr <- q3 - q1
      if (iqr > 0) outlier_flags <- outlier_flags + as.numeric(abs(xj - median(xj, na.rm = TRUE)) > 3 * iqr)
    }
    cov_outlier_rate <<- mean(outlier_flags > 0)
  }, error = function(e) { cov_outlier_rate <<- 0 })
  
  y_mad <- mad(data$Y); y_median <- median(data$Y)
  y_outlier_rate <- ifelse(y_mad > 0, mean(abs(data$Y - y_median) > 3 * y_mad), 0)
  
  ps_formula <- as.formula(paste("T ~", paste(names(data)[-(1:2)], collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial())
  ps_hat <- predict(ps_model, type = "response")
  ps_overlap <- mean(ps_hat > 0.1 & ps_hat < 0.9)
  
  treat_imbalance <- abs(mean(data$T) - 0.5)
  w1 <- data$T / pmax(0.01, ps_hat); w0 <- (1 - data$T) / pmax(0.01, 1 - ps_hat)
  weight_cv <- sd(c(w1, w0)) / mean(c(w1, w0))
  
  y_kurt <- 0
  tryCatch({
    if (require("moments", quietly = TRUE)) {
      y_kurt <- moments::kurtosis(data$Y) - 3
    } else {
      y_scaled <- (data$Y - mean(data$Y)) / sd(data$Y); y_kurt <- mean(y_scaled^4) - 3
    }
  }, error = function(e) { y_kurt <<- 0 })
  
  contamination_score <- cov_outlier_rate + y_outlier_rate + (1 - ps_overlap) +
    treat_imbalance + min(weight_cv / 10, 1) + max(0, y_kurt / 10)
  
  if (contamination_score < 0.5) {
    selected <- "reg_adj"; reason <- "Low contamination"; confidence <- "high"
  } else if (cov_outlier_rate > 0.2 & y_outlier_rate > 0.2) {
    selected <- "aipw_robust"; reason <- "Severe mixed contamination"; confidence <- "high"
  } else if (cov_outlier_rate > y_outlier_rate & cov_outlier_rate > 0.1) {
    selected <- "dr_robust"; reason <- "Covariate outliers dominant"; confidence <- "medium"
  } else if (y_outlier_rate > 0.15) {
    selected <- "aipw_robust"; reason <- "Severe outcome contamination"; confidence <- "high"
  } else if (weight_cv > 5) {
    selected <- "ipw_robust"; reason <- "Extreme weights"; confidence <- "medium"
  } else {
    selected <- "ipw_robust"; reason <- "Moderate contamination"; confidence <- "medium"
  }
  
  list(selected = selected, reason = reason, confidence = confidence,
       diagnostics = c(cov_outlier_rate = cov_outlier_rate, y_outlier_rate = y_outlier_rate,
                       ps_overlap = ps_overlap, treat_imbalance = treat_imbalance,
                       weight_cv = weight_cv, y_kurt = y_kurt, contamination_score = contamination_score))
}

# =============================================================================
# 6. REAL DATA: LALONDE (MatchIt uyumlu, degismedi)
# =============================================================================
load_lalonde_data <- function() {
  cat("\n--- LALONDE DATASET ---\n")
  lalonde_loaded <- FALSE
  
  if (require("MatchIt", quietly = TRUE)) {
    tryCatch({
      lalonde <- MatchIt::lalonde
      if (nrow(lalonde) > 0) { cat("Loaded from MatchIt::lalonde\n"); lalonde_loaded <- TRUE }
    }, error = function(e) { cat("MatchIt found but data load failed.\n") })
  }
  if (!lalonde_loaded && require("Matching", quietly = TRUE)) {
    tryCatch({ data(lalonde); if (exists("lalonde") && nrow(lalonde) > 0) lalonde_loaded <- TRUE }, error = function(e) {})
  }
  if (!lalonde_loaded) {
    cat("Creating built-in LaLonde data...\n")
    set.seed(123)
    n_treat <- 185; n_ctrl <- 429
    treat_data <- data.frame(age = round(rnorm(n_treat, 25.8, 7.2)), educ = round(rnorm(n_treat, 10.3, 2.0)),
                             black = rbinom(n_treat, 1, 0.84), hispan = rbinom(n_treat, 1, 0.06), married = rbinom(n_treat, 1, 0.19),
                             nodegree = rbinom(n_treat, 1, 0.71), re74 = rlnorm(n_treat, 7.2, 1.8), re75 = rlnorm(n_treat, 7.5, 1.6),
                             re78 = rlnorm(n_treat, 8.6, 1.4), treat = 1)
    ctrl_data <- data.frame(age = round(rnorm(n_ctrl, 34.9, 10.4)), educ = round(rnorm(n_ctrl, 12.1, 3.1)),
                            black = rbinom(n_ctrl, 1, 0.45), hispan = rbinom(n_ctrl, 1, 0.12), married = rbinom(n_ctrl, 1, 0.51),
                            nodegree = rbinom(n_ctrl, 1, 0.31), re74 = rlnorm(n_ctrl, 8.4, 1.2), re75 = rlnorm(n_ctrl, 8.5, 1.1),
                            re78 = rlnorm(n_ctrl, 8.8, 1.2), treat = 0)
    lalonde <- rbind(treat_data, ctrl_data)
  }
  
  col_names <- names(lalonde)
  if ("race" %in% col_names) {
    black <- as.numeric(lalonde$race == "black"); hispan <- as.numeric(lalonde$race == "hispan")
  } else {
    black <- if ("black" %in% col_names) lalonde$black else 0
    hispan <- if ("hispan" %in% col_names) lalonde$hispan else if ("hisp" %in% col_names) lalonde$hisp else 0
  }
  
  data_prep <- data.frame(Y = lalonde$re78, T = lalonde$treat, age = lalonde$age, educ = lalonde$educ,
                          black = black, hispan = hispan, married = lalonde$married, nodegree = lalonde$nodegree,
                          re74 = lalonde$re74, re75 = lalonde$re75)
  data_prep <- na.omit(data_prep)
  cat(sprintf("Final n=%d, Treatment=%d (%.1f%%)\n", nrow(data_prep), sum(data_prep$T), 100*mean(data_prep$T)))
  data_prep
}

analyze_lalonde <- function() {
  data <- load_lalonde_data()
  est_naive <- tryCatch(calc_ipw_naive(data), error = function(e) NA)
  est_ipw_r <- tryCatch(calc_ipw_robust(data), error = function(e) NA)
  est_dr_r  <- tryCatch(calc_dr_robust(data), error = function(e) NA)
  est_aipw_r<- tryCatch(calc_aipw_robust(data), error = function(e) NA)
  est_reg   <- tryCatch(calc_reg_adj(data), error = function(e) NA)
  adaptive <- tryCatch(select_estimator_adaptive(data), error = function(e) list(selected="error", reason="error", confidence="low"))
  est_adaptive <- switch(adaptive$selected, "ipw_robust" = est_ipw_r, "dr_robust" = est_dr_r,
                         "aipw_robust" = est_aipw_r, "reg_adj" = est_reg, est_ipw_r)
  results <- data.frame(Method = c("Naive IPW","Robust IPW","DR Robust","AIPW Robust","Reg Adj","Adaptive"),
                        ATE_Estimate = c(est_naive, est_ipw_r, est_dr_r, est_aipw_r, est_reg, est_adaptive))
  cat("\nLaLonde Results:\n"); print(results)
  cat(sprintf("Adaptive: %s | %s\n", adaptive$selected, adaptive$reason))
  list(results = results, adaptive = adaptive, data = data)
}

# =============================================================================
# 7. REAL DATA: NHANES (synthetic, degismedi)
# =============================================================================
load_nhanes_data <- function() {
  cat("\n--- NHANES (Synthetic) ---\n")
  set.seed(456); n <- 5000
  d <- data.frame(age = round(rnorm(n, 45, 15)), gender = rbinom(n, 1, 0.5),
                  race = sample(1:4, n, prob = c(0.7,0.12,0.08,0.10), replace = TRUE), bmi = rnorm(n, 28, 6),
                  systolic = rnorm(n, 120, 15), diastolic = rnorm(n, 78, 10), cholesterol = rnorm(n, 200, 35),
                  diabetes = rbinom(n, 1, 0.12), smoking = rbinom(n, 1, 0.22))
  d$health <- 100 - 0.3*d$age - 2*d$bmi - 0.5*d$systolic - 5*d$diabetes - 8*d$smoking + rnorm(n, 0, 10)
  dp <- data.frame(Y = d$health, T = d$smoking, age = d$age, gender = d$gender, race = d$race,
                   bmi = d$bmi, systolic = d$systolic, diastolic = d$diastolic,
                   cholesterol = d$cholesterol, diabetes = d$diabetes)
  na.omit(dp)
}

analyze_nhanes <- function() {
  data <- load_nhanes_data()
  est_naive <- tryCatch(calc_ipw_naive(data), error = function(e) NA)
  est_ipw_r <- tryCatch(calc_ipw_robust(data), error = function(e) NA)
  est_dr_r  <- tryCatch(calc_dr_robust(data), error = function(e) NA)
  est_aipw_r<- tryCatch(calc_aipw_robust(data), error = function(e) NA)
  est_reg   <- tryCatch(calc_reg_adj(data), error = function(e) NA)
  adaptive <- tryCatch(select_estimator_adaptive(data), error = function(e) list(selected="error", reason="error", confidence="low"))
  est_adaptive <- switch(adaptive$selected, "ipw_robust" = est_ipw_r, "dr_robust" = est_dr_r,
                         "aipw_robust" = est_aipw_r, "reg_adj" = est_reg, est_ipw_r)
  results <- data.frame(Method = c("Naive IPW","Robust IPW","DR Robust","AIPW Robust","Reg Adj","Adaptive"),
                        ATE_Estimate = c(est_naive, est_ipw_r, est_dr_r, est_aipw_r, est_reg, est_adaptive))
  cat("\nNHANES Results:\n"); print(results)
  list(results = results, adaptive = adaptive, data = data)
}

# =============================================================================
# 8. MC ENGINE
# =============================================================================
run_single_mc <- function(n, p, true_ate, contamination_type, contamination_rate, contamination_strength) {
  data <- generate_data(n, p, true_ate, contamination_type, contamination_rate, contamination_strength)$data
  est_naive<-NA; est_ipw_robust<-NA; est_dr_robust<-NA; est_aipw_robust<-NA; est_reg_adj<-NA; est_adaptive<-NA
  adaptive_selected <- "error"; adaptive_confidence <- "low"
  tryCatch({ est_naive <- calc_ipw_naive(data) }, error = function(e) {})
  tryCatch({ est_ipw_robust <- calc_ipw_robust(data) }, error = function(e) {})
  tryCatch({ est_dr_robust <- calc_dr_robust(data) }, error = function(e) {})
  tryCatch({ est_aipw_robust <- calc_aipw_robust(data) }, error = function(e) {})
  tryCatch({ est_reg_adj <- calc_reg_adj(data) }, error = function(e) {})
  tryCatch({
    adaptive <- select_estimator_adaptive(data)
    adaptive_selected <- adaptive$selected; adaptive_confidence <- adaptive$confidence
    est_adaptive <- switch(adaptive_selected, "ipw_robust" = est_ipw_robust, "dr_robust" = est_dr_robust,
                           "aipw_robust" = est_aipw_robust, "reg_adj" = est_reg_adj, est_ipw_robust)
  }, error = function(e) {})
  c(naive = est_naive, ipw_robust = est_ipw_robust, dr_robust = est_dr_robust,
    aipw_robust = est_aipw_robust, reg_adj = est_reg_adj, adaptive = est_adaptive,
    selected = adaptive_selected, confidence = adaptive_confidence)
}

estimate_time <- function(elapsed, completed, total) {
  if (completed == 0) return("Unknown")
  remaining <- (elapsed / completed) * (total - completed)
  sprintf("%d saat %d dakika", floor(remaining/3600), floor((remaining %% 3600)/60))
}

# =============================================================================
# 9. PRODUCTION MC RUN - [FIX 3] header bug giderildi, [FIX 4] n_mc kontrolu
# =============================================================================
run_monte_carlo_parallel <- function(config = CONFIG) {
  
  n_cores <- min(config$n_cores, detectCores() - 1)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  if (USE_DORNG) registerDoRNG(config$seed)   # [FIX 5] bagimsiz/tekrarlanabilir RNG akislari
  cat(sprintf("\n%d cekirdek kullaniliyor%s\n", n_cores, ifelse(USE_DORNG, " (doRNG aktif)", " (UYARI: doRNG YOK)")))
  
  clusterExport(cl, c("generate_data", "calc_ipw_naive", "calc_ipw_robust", "calc_dr_robust",
                      "calc_aipw_robust", "calc_reg_adj", "select_estimator_adaptive",
                      "run_single_mc", "BETA_PS_FULL", "BETA_OUT_FULL"))
  clusterEvalQ(cl, { library(MASS); library(robustbase) })
  
  scenarios <- expand.grid(n = config$n_values, p = config$p_values,
                           contamination_type = config$contamination_types,
                           contamination_rate = config$contamination_rates,
                           contamination_strength = config$contamination_strengths,
                           stringsAsFactors = FALSE)
  total_scenarios <- nrow(scenarios)
  cat(sprintf("Toplam senaryo: %d | Tekrar/senaryo: %d | Toplam simulasyon: %d\n",
              total_scenarios, config$n_mc, total_scenarios * config$n_mc))
  
  checkpoint_file <- config$checkpoint_file
  completed_scenario_ids <- integer(0)
  
  if (file.exists(checkpoint_file)) {
    existing <- tryCatch(read.csv(checkpoint_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(existing) && nrow(existing) > 0) {
      # [FIX 4] KRITIK GUVENLIK KONTROLU: checkpoint farkli n_mc ile mi olusturulmus?
      existing_nmc <- unique(existing$n_mc)
      if (!all(existing_nmc == config$n_mc)) {
        stop(sprintf(paste0("DURDURULDU: '%s' dosyasindaki checkpoint n_mc=%s ile olusturulmus, ",
                            "ama CONFIG$n_mc=%d. Bu dosyayi SILIN veya YENIDEN ADLANDIRIN, ",
                            "yoksa dusuk replikasyonlu eski sonuclar 'tamamlandi' sanilip atlanir."),
                     checkpoint_file, paste(existing_nmc, collapse=","), config$n_mc))
      }
      completed_scenario_ids <- unique(existing$scenario_id)
      cat(sprintf("Checkpoint bulundu: %d senaryo tamamlanmis (n_mc dogrulandi). Devam ediliyor...\n",
                  length(completed_scenario_ids)))
    }
  }
  
  results <- list()
  start_time <- Sys.time()
  n_done_this_session <- 0
  
  for (i in 1:total_scenarios) {
    if (i %in% completed_scenario_ids) next
    
    scen <- scenarios[i, ]
    if (config$verbose) {
      cat(sprintf("\n[%d/%d] n=%d p=%d %s rate=%.2f strength=%.1f\n",
                  i, total_scenarios, scen$n, scen$p, scen$contamination_type,
                  scen$contamination_rate, scen$contamination_strength))
    }
    
    scen_start <- Sys.time()
    if (USE_DORNG) {
      mc_results <- foreach(rep = 1:config$n_mc, .combine = rbind,
                            .packages = c("MASS", "robustbase")) %dorng% {
                              res <- run_single_mc(scen$n, scen$p, config$true_ate, scen$contamination_type,
                                                   scen$contamination_rate, scen$contamination_strength)
                              c(as.numeric(res[1:6]), as.character(res[7]), as.character(res[8]))
                            }
    } else {
      mc_results <- foreach(rep = 1:config$n_mc, .combine = rbind,
                            .packages = c("MASS", "robustbase")) %dopar% {
                              res <- run_single_mc(scen$n, scen$p, config$true_ate, scen$contamination_type,
                                                   scen$contamination_rate, scen$contamination_strength)
                              c(as.numeric(res[1:6]), as.character(res[7]), as.character(res[8]))
                            }
    }
    scen_elapsed <- as.numeric(difftime(Sys.time(), scen_start, units = "secs"))
    
    methods <- c("naive", "ipw_robust", "dr_robust", "aipw_robust", "reg_adj", "adaptive")
    for (j in seq_along(methods)) {
      estimates <- as.numeric(mc_results[, j]); estimates <- estimates[!is.na(estimates)]
      if (length(estimates) > 0) {
        bias <- mean(estimates) - config$true_ate; mse <- mean((estimates - config$true_ate)^2)
        result_row <- data.frame(scenario_id = i, method = methods[j], n = scen$n, p = scen$p,
                                 contamination_type = scen$contamination_type, contamination_rate = scen$contamination_rate,
                                 contamination_strength = scen$contamination_strength, true_ate = config$true_ate, n_mc = config$n_mc,
                                 bias = bias, mse = mse, var = var(estimates), rmse = sqrt(mse), mae = mean(abs(estimates - config$true_ate)),
                                 coverage = mean(abs(estimates - config$true_ate) < 1.96 * sd(estimates)),
                                 n_valid = length(estimates), elapsed_sec = scen_elapsed, stringsAsFactors = FALSE)
      } else {
        result_row <- data.frame(scenario_id = i, method = methods[j], n = scen$n, p = scen$p,
                                 contamination_type = scen$contamination_type, contamination_rate = scen$contamination_rate,
                                 contamination_strength = scen$contamination_strength, true_ate = config$true_ate, n_mc = config$n_mc,
                                 bias = NA, mse = NA, var = NA, rmse = NA, mae = NA, coverage = NA,
                                 n_valid = 0, elapsed_sec = scen_elapsed, stringsAsFactors = FALSE)
        cat(sprintf("  UYARI: %s icin gecerli tahmin yok\n", methods[j]))
      }
      results[[length(results) + 1]] <- result_row
    }
    
    # [FIX 3] Basitlestirilmis, guvenli checkpoint yazimi - stub dosya YOK
    if (length(results) > 0) {
      batch <- do.call(rbind, results)
      write.table(batch, checkpoint_file, append = file.exists(checkpoint_file),
                  sep = ",", row.names = FALSE, col.names = !file.exists(checkpoint_file))
      results <- list()
    }
    
    n_done_this_session <- n_done_this_session + 1
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    remaining <- estimate_time(elapsed, n_done_this_session,
                               total_scenarios - length(completed_scenario_ids))
    if (config$verbose) {
      cat(sprintf("  Senaryo suresi: %.1f sn | Bu oturumda gecen: %.1f dk | Tahmini kalan: %s\n",
                  scen_elapsed, elapsed/60, remaining))
    }
  }
  
  stopCluster(cl)
  cat("\n--- SONUCLAR DERLENIYOR ---\n")
  final_results <- read.csv(checkpoint_file, stringsAsFactors = FALSE)
  write.csv(final_results, config$results_file, row.names = FALSE)
  cat(sprintf("Kaydedildi: %s (%d satir)\n", config$results_file, nrow(final_results)))
  final_results
}

# =============================================================================
# 10. VISUALIZATION (degismedi)
# =============================================================================
plot_performance <- function(results) {
  p1 <- ggplot(results, aes(x = factor(contamination_rate), y = mse, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(contamination_type ~ contamination_strength + n, scales = "free_y") +
    labs(title = "MSE Comparison", x = "Contamination Rate", y = "MSE") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- ggplot(results, aes(x = factor(contamination_rate), y = abs(bias), fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(contamination_type ~ contamination_strength + n, scales = "free_y") +
    labs(title = "Absolute Bias", x = "Contamination Rate", y = "|Bias|") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p3 <- ggplot(results, aes(x = factor(contamination_rate), y = rmse, color = method, group = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(contamination_type ~ contamination_strength + n, scales = "free_y") +
    labs(title = "RMSE", x = "Contamination Rate", y = "RMSE") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p4 <- ggplot(results, aes(x = factor(n), y = mse, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") + facet_grid(contamination_type ~ p, scales = "free_y") +
    labs(title = "MSE by n", x = "n", y = "MSE") + theme_minimal()
  list(mse = p1, bias = p2, rmse = p3, sample_size = p4)
}

plot_adaptive_analysis <- function(results) {
  ar <- results[results$method == "adaptive", ]
  p1 <- ggplot(ar, aes(x = factor(contamination_rate), y = rmse, color = factor(contamination_strength))) +
    geom_line(linewidth = 1) + geom_point(size = 3) + facet_grid(contamination_type ~ n, scales = "free_y") +
    labs(title = "Adaptive RMSE", x = "Rate", y = "RMSE", color = "Strength") + theme_minimal()
  mc <- results[results$method %in% c("naive", "adaptive"), ]
  p2 <- ggplot(mc, aes(x = factor(contamination_rate), y = rmse, color = method, group = method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2.5) + facet_grid(contamination_type ~ n, scales = "free_y") +
    labs(title = "Adaptive vs Naive", x = "Rate", y = "RMSE") + theme_minimal()
  list(adaptive_rmse = p1, comparison = p2)
}

# =============================================================================
# 11. MAIN EXECUTION
# =============================================================================
cat(rep("=", 70), "\n", sep = "")
cat("ADAPTIVE MONTE CARLO FRAMEWORK - PRODUCTION (n_mc=2000)\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("STEP 1: THEORETICAL ANALYSIS\n")
theory_data <- generate_data(n = 500, p = 5, true_ate = 2.0, contamination_type = "mixed",
                             contamination_rate = 0.10, contamination_strength = 5)
if_result <- plot_influence_functions(theory_data$data)
ggsave("theory_if_plot.pdf", if_result$plot, width = 10, height = 6)
bias_result <- calc_asymptotic_bias(generate_data(500, 5, 2.0, "none", 0, 0)$data, theory_data$data, epsilon = 0.10)
bp_result <- compare_breakdown_points(theory_data$data)
ggsave("theory_breakdown_plot.pdf", bp_result$plot, width = 10, height = 6)

cat("\n\nSTEP 2: REAL DATA APPLICATIONS\n")
lalonde_result <- analyze_lalonde()
nhanes_result <- analyze_nhanes()

cat("\n\nSTEP 3: MONTE CARLO SIMULATION (n_mc=2000, TAM KAPSAM)\n")
cat(sprintf("Tahmini toplam sure: ~9.9 gun (675 senaryo, gercek n_mc=10 olcumune dayali)\n"))
mc_results <- run_monte_carlo_parallel(CONFIG)

cat("\n\nSTEP 4: VISUALIZATION AND SUMMARY\n")
plots <- plot_performance(mc_results)
adaptive_plots <- plot_adaptive_analysis(mc_results)
pdf("mc_full_results_production.pdf", width = 14, height = 10)
print(plots$mse); print(plots$bias); print(plots$rmse); print(plots$sample_size)
print(adaptive_plots$adaptive_rmse); print(adaptive_plots$comparison)
dev.off()

summary_table <- aggregate(cbind(bias, mse, rmse, mae) ~ method + contamination_type, data = mc_results, FUN = mean)
cat("\n--- OVERALL PERFORMANCE SUMMARY ---\n")
print(summary_table)

cat("\n", rep("=", 70), "\n", sep = "")
cat("PRODUCTION RUN TAMAMLANDI\n")
cat(rep("=", 70), "\n", sep = "")

# =============================================================================
# CALISTIRMA TALIMATLARI (script disinda, sadece bilgi amacli)
# =============================================================================
# Windows'ta komut satirindan (RStudio DEGIL) calistirmak icin:
#   1) Komut Istemi / PowerShell acin, script'in oldugu klasore gidin
#   2) "C:\Program Files\R\R-4.x.x\bin\Rscript.exe" production_run.R > production_log.txt 2>&1
#   3) Guc Secenekleri > Uyku/Ekran kapanmasi -> "Asla" yapin
#   4) Ilerlemeyi kontrol etmek icin production_log.txt dosyasini herhangi
#      bir zamanda acip son satirlara bakabilirsiniz (dosya surekli guncellenir)
#   5) Kesinti olursa (elektrik, yeniden baslatma vb.) AYNI KOMUTU tekrar
#      calistirmaniz yeterli - checkpoint kaldigi yerden devam eder
#      (test edildi, dogrulandi).