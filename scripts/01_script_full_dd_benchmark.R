# ============================================================
# FULL MASTER DD BENCHMARK SCRIPT
# Robust reconstruction of the original benchmark macro
# Ready for Posit / RStudio
# Date: 2026-03-19
# ============================================================

rm(list = ls())

pkg <- c("ddalpha", "MASS", "e1071", "RSSL", "dplyr", "kernlab")
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
invisible(lapply(pkg, require, character.only = TRUE))

set.seed(42)

OUTPUT_DIR <- "outputs"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

DATASETS <- c(
  "baby", "banknoten", "biomed", "bloodtransfusion",
  "breast_cancer_wisconsin", "bupa",
  "chemdiab_1vs2", "chemdiab_1vs3", "chemdiab_2vs3",
  "cloud", "crabB_MvsF", "crabF_BvsO", "crabM_BvsO", "crabO_MvsF",
  "crab_BvsO", "crab_MvsF", "cricket_CvsP",
  "ecoli_cpvsim", "ecoli_cpvspp", "ecoli_imvspp",
  "gemsen_MvsF", "glass", "groessen_MvsF", "haberman", "heart",
  "hemophilia", "indian_liver_patient_1vs2", "indian_liver_patient_FvsM",
  "iris_setosavsversicolor", "iris_setosavsvirginica", "iris_versicolorvsvirginica",
  "irish_ed_MvsF", "kidney", "pima", "plasma_retinol_MvsF",
  "segmentation", "socmob_IvsNI", "socmob_WvsB", "tae", "tennis_MvsF",
  "tips_DvsN", "tips_MvsF", "uscrime_SvsN", "vertebral_column",
  "veteran_lung_cancer", "vowel_MvsF", "wine_1vs2", "wine_1vs3", "wine_2vs3"
)

UNSUP_PROP <- 0.15
USE_SEMI_SUPERVISED <- TRUE
DEPTH_METHOD <- "halfspace"
VERBOSE <- TRUE
FAST_MODE <- F
MAX_FOLDS <- 2
QUIET_DDALPHA <- TRUE

log_msg <- function(...) {
  if (isTRUE(VERBOSE)) cat(..., "\n")
}

safe_factor01 <- function(y) {
  y <- as.factor(y)
  if (length(levels(y)) != 2) stop("Response is not binary.")
  levels(y) <- c("0", "1")
  y
}

prepare_dataset <- function(dataset_name) {
  data <- getdata(dataset_name)
  colnames(data)[ncol(data)] <- "C"
  
  if (dataset_name %in% c("cloud", "gemsen_MvsF", "wine_2vs3")) {
    data$C <- as.numeric(data$C)
  }
  
  data <- na.omit(data)
  X <- subset(data, select = -C)
  X <- as.data.frame(lapply(X, as.numeric))
  y <- safe_factor01(data$C)
  
  if (!all(sapply(X, is.numeric))) stop("Not all predictors are numeric.")
  if (nrow(X) < 5) stop("Too few observations.")
  if (any(table(y) < 2)) stop("A class has fewer than 2 observations.")
  
  list(X = X, y = y)
}

scale_train_test <- function(X_train, X_test) {
  mu <- sapply(X_train, mean)
  sdv <- sapply(X_train, sd)
  sdv[is.na(sdv) | sdv == 0] <- 1
  list(
    train = as.data.frame(scale(X_train, center = mu, scale = sdv)),
    test  = as.data.frame(scale(X_test, center = mu, scale = sdv))
  )
}

compute_dd_representation <- function(X_train, y_train, X_test,
                                      depth_method = c("halfspace", "Mahalanobis")) {
  depth_method <- match.arg(depth_method)
  
  X_train_m <- as.matrix(X_train)
  X_test_m  <- as.matrix(X_test)
  
  ref0 <- X_train_m[y_train == "0", , drop = FALSE]
  ref1 <- X_train_m[y_train == "1", , drop = FALSE]
  
  if (nrow(ref0) < 2 || nrow(ref1) < 2) {
    stop("One class has fewer than 2 training observations.")
  }
  
  depth_fun <- switch(
    depth_method,
    "halfspace" = function(x, data) ddalpha::depth.halfspace(x = x, data = data),
    "Mahalanobis" = function(x, data) ddalpha::depth.Mahalanobis(x = x, data = data)
  )
  
  D1_train <- sapply(seq_len(nrow(X_train_m)), function(i)
    depth_fun(X_train_m[i, , drop = FALSE], ref0)
  )
  D2_train <- sapply(seq_len(nrow(X_train_m)), function(i)
    depth_fun(X_train_m[i, , drop = FALSE], ref1)
  )
  D1_test <- depth_fun(X_test_m, ref0)
  D2_test <- depth_fun(X_test_m, ref1)
  
  dd_train <- data.frame(D1 = D1_train, D2 = D2_train)
  dd_test  <- data.frame(D1 = D1_test, D2 = D2_test)
  
  ddpoly_train <- transform(dd_train, D1_sq = D1^2, D2_sq = D2^2, D1D2 = D1 * D2)
  ddpoly_test  <- transform(dd_test,  D1_sq = D1^2, D2_sq = D2^2, D1D2 = D1 * D2)
  
  list(dd_train = dd_train, dd_test = dd_test,
       ddpoly_train = ddpoly_train, ddpoly_test = ddpoly_test)
}

drop_constant_cols <- function(df) {
  keep <- sapply(df, function(x) length(unique(x)) > 1)
  df[, keep, drop = FALSE]
}

drop_lda_bad_cols <- function(df, y) {
  keep <- sapply(seq_along(df), function(j) {
    x <- df[[j]]
    vars_by_class <- tapply(x, y, function(v) var(v, na.rm = TRUE))
    all(!is.na(vars_by_class) & vars_by_class > 0)
  })
  df[, keep, drop = FALSE]
}

fit_predict_supervised <- function(method, X_train, y_train, X_test, dd_obj) {
  if (method == "LDA.orig") {
    fit <- MASS::lda(x = X_train, grouping = y_train)
    return(as.character(predict(fit, X_test)$class)[1])
  }
  
  if (method == "SVM.orig") {
    fit <- e1071::svm(x = X_train, y = y_train, kernel = "radial", scale = FALSE)
    return(as.character(predict(fit, X_test))[1])
  }
  
  if (method == "DD.alpha") {
    assign(".ddalpha_tmp_train", data.frame(C = y_train, X_train, check.names = FALSE), envir = .GlobalEnv)
    assign(".ddalpha_tmp_test",  data.frame(X_test, check.names = FALSE), envir = .GlobalEnv)
    on.exit({
      if (exists(".ddalpha_tmp_train", envir = .GlobalEnv)) rm(".ddalpha_tmp_train", envir = .GlobalEnv)
      if (exists(".ddalpha_tmp_test",  envir = .GlobalEnv)) rm(".ddalpha_tmp_test",  envir = .GlobalEnv)
    }, add = TRUE)
    
    if (QUIET_DDALPHA) {
      fit <- suppressMessages(suppressWarnings(
        ddalpha::ddalpha.train(C ~ ., data = .ddalpha_tmp_train,
                               depth = "halfspace", separator = "alpha")
      ))
    } else {
      fit <- ddalpha::ddalpha.train(C ~ ., data = .ddalpha_tmp_train,
                                    depth = "halfspace", separator = "alpha")
    }
    return(as.character(predict(fit, .ddalpha_tmp_test))[1])
  }
  
  if (method == "LDA.DD") {
    tr <- drop_lda_bad_cols(dd_obj$dd_train, y_train)
    if (ncol(tr) == 0) stop("No valid DD columns for LDA after removing within-class constants.")
    te <- dd_obj$dd_test[, colnames(tr), drop = FALSE]
    fit <- MASS::lda(x = tr, grouping = y_train)
    return(as.character(predict(fit, te)$class)[1])
  }
  
  if (method == "SVM.DD") {
    fit <- e1071::svm(x = dd_obj$dd_train, y = y_train, kernel = "radial", scale = FALSE)
    return(as.character(predict(fit, dd_obj$dd_test))[1])
  }
  
  if (method == "LDA.poly") {
    tr <- drop_lda_bad_cols(dd_obj$ddpoly_train, y_train)
    if (ncol(tr) == 0) stop("No valid DD-poly columns for LDA after removing within-class constants.")
    te <- dd_obj$ddpoly_test[, colnames(tr), drop = FALSE]
    fit <- MASS::lda(x = tr, grouping = y_train)
    return(as.character(predict(fit, te)$class)[1])
  }
  
  if (method == "SVM.poly") {
    fit <- e1071::svm(x = dd_obj$ddpoly_train, y = y_train, kernel = "radial", scale = FALSE)
    return(as.character(predict(fit, dd_obj$ddpoly_test))[1])
  }
  
  stop("Unknown supervised method")
}

mask_labels <- function(y, prop = 0.15) {
  y_ssl <- as.character(y)
  n_mask <- max(1, floor(length(y_ssl) * prop))
  idx_mask <- sample(seq_along(y_ssl), size = n_mask, replace = FALSE)
  y_ssl[idx_mask] <- NA
  factor(y_ssl, levels = c("0", "1"))
}

build_ssl_problem <- function(X_train, y_train_ssl) {
  df_ssl <- data.frame(C = y_train_ssl, X_train)
  RSSL::df_to_matrices(df_ssl, C ~ .)
}

fit_predict_ssl <- function(method, X_train_sc, X_test_sc, y_train_ssl, dd_obj) {
  tab_labeled <- table(na.omit(y_train_ssl))
  if (length(tab_labeled) < 2 || any(tab_labeled < 2)) {
    stop("Too few labeled observations per class after masking.")
  }
  
  if (method == "SSL.orig") {
    problem <- build_ssl_problem(X_train_sc, y_train_ssl)
    fit <- RSSL::LaplacianKernelLeastSquaresClassifier(
      X = problem$X, y = problem$y, X_u = problem$X_u,
      kernel = kernlab::rbfdot(0.5), lambda = 0.001, gamma = 100,
      normalized_laplacian = TRUE, scale = TRUE
    )
    return(as.character(predict(fit, as.matrix(X_test_sc)))[1])
  }
  
  if (method == "SSL.DD") {
    problem <- build_ssl_problem(dd_obj$dd_train, y_train_ssl)
    fit <- RSSL::LaplacianKernelLeastSquaresClassifier(
      X = problem$X, y = problem$y, X_u = problem$X_u,
      kernel = kernlab::rbfdot(0.5), lambda = 0.001, gamma = 100,
      normalized_laplacian = TRUE, scale = TRUE
    )
    return(as.character(predict(fit, as.matrix(dd_obj$dd_test)))[1])
  }
  
  if (method == "SSL.poly") {
    tr <- drop_constant_cols(dd_obj$ddpoly_train)
    if (ncol(tr) == 0) stop("All DD-poly columns constant.")
    te <- dd_obj$ddpoly_test[, colnames(tr), drop = FALSE]
    problem <- build_ssl_problem(tr, y_train_ssl)
    fit <- RSSL::LaplacianKernelLeastSquaresClassifier(
      X = problem$X, y = problem$y, X_u = problem$X_u,
      kernel = kernlab::rbfdot(0.5), lambda = 0.001, gamma = 100,
      normalized_laplacian = TRUE, scale = TRUE
    )
    return(as.character(predict(fit, as.matrix(te)))[1])
  }
  
  stop("Unknown SSL method")
}

evaluate_dataset <- function(dataset_name) {
  methods_super <- c("LDA.orig", "SVM.orig", "DD.alpha", "LDA.DD", "SVM.DD", "LDA.poly", "SVM.poly")
  methods_ssl <- c("SSL.orig", "SSL.DD", "SSL.poly")
  
  prepared <- prepare_dataset(dataset_name)
  X <- prepared$X
  y <- prepared$y
  n <- nrow(X)
  
  fold_seq <- if (FAST_MODE) seq_len(min(n, MAX_FOLDS)) else seq_len(n)
  
  super_attempted <- setNames(rep(0L, length(methods_super)), methods_super)
  super_successful <- setNames(rep(0L, length(methods_super)), methods_super)
  super_correct <- setNames(rep(0L, length(methods_super)), methods_super)
  super_first_error <- setNames(as.list(rep(NA_character_, length(methods_super))), methods_super)
  
  ssl_attempted <- setNames(rep(0L, length(methods_ssl)), methods_ssl)
  ssl_successful <- setNames(rep(0L, length(methods_ssl)), methods_ssl)
  ssl_correct <- setNames(rep(0L, length(methods_ssl)), methods_ssl)
  ssl_first_error <- setNames(as.list(rep(NA_character_, length(methods_ssl))), methods_ssl)
  
  for (i in fold_seq) {
    if (VERBOSE && (i == 1 || i %% 50 == 0)) {
      log_msg(sprintf("Dataset %s - Fold %d / %d", dataset_name, i, max(fold_seq)))
    }
    
    X_train <- X[-i, , drop = FALSE]
    X_test  <- X[i, , drop = FALSE]
    y_train <- y[-i]
    y_test  <- y[i]
    
    sc <- scale_train_test(X_train, X_test)
    X_train_sc <- sc$train
    X_test_sc  <- sc$test
    dd_obj <- compute_dd_representation(X_train_sc, y_train, X_test_sc, depth_method = DEPTH_METHOD)
    
    for (m in methods_super) {
      super_attempted[m] <- super_attempted[m] + 1L
      res <- tryCatch({
        pred <- fit_predict_supervised(m, X_train_sc, y_train, X_test_sc, dd_obj)
        list(ok = TRUE, pred = pred, err = NA_character_)
      }, error = function(e) {
        list(ok = FALSE, pred = NA_character_, err = conditionMessage(e))
      })
      
      if (isTRUE(res$ok) && !is.na(res$pred)) {
        super_successful[m] <- super_successful[m] + 1L
        super_correct[m] <- super_correct[m] + as.integer(res$pred == as.character(y_test))
      } else {
        if (is.na(super_first_error[[m]])) super_first_error[[m]] <- res$err
      }
    }
    
    if (USE_SEMI_SUPERVISED) {
      y_train_ssl <- mask_labels(y_train, prop = UNSUP_PROP)
      
      for (m in methods_ssl) {
        ssl_attempted[m] <- ssl_attempted[m] + 1L
        res <- tryCatch({
          pred <- fit_predict_ssl(m, X_train_sc, X_test_sc, y_train_ssl, dd_obj)
          list(ok = TRUE, pred = pred, err = NA_character_)
        }, error = function(e) {
          list(ok = FALSE, pred = NA_character_, err = conditionMessage(e))
        })
        
        if (isTRUE(res$ok) && !is.na(res$pred)) {
          ssl_successful[m] <- ssl_successful[m] + 1L
          ssl_correct[m] <- ssl_correct[m] + as.integer(res$pred == as.character(y_test))
        } else {
          if (is.na(ssl_first_error[[m]])) ssl_first_error[[m]] <- res$err
        }
      }
    }
  }
  
  out <- list()
  for (m in methods_super) {
    out[[length(out) + 1]] <- data.frame(
      dataset = dataset_name,
      method = m,
      setting = "supervised",
      attempted_folds = super_attempted[m],
      successful_folds = super_successful[m],
      failed_folds = super_attempted[m] - super_successful[m],
      accuracy = if (super_successful[m] > 0) super_correct[m] / super_successful[m] else NA_real_,
      correct = super_correct[m],
      status = if (super_successful[m] == super_attempted[m] && super_attempted[m] > 0) "OK" else if (super_successful[m] > 0) "PARTIAL" else "FAIL",
      first_error_message = as.character(super_first_error[[m]]),
      stringsAsFactors = FALSE
    )
  }
  
  if (USE_SEMI_SUPERVISED) {
    for (m in methods_ssl) {
      out[[length(out) + 1]] <- data.frame(
        dataset = dataset_name,
        method = m,
        setting = "semi_supervised",
        attempted_folds = ssl_attempted[m],
        successful_folds = ssl_successful[m],
        failed_folds = ssl_attempted[m] - ssl_successful[m],
        accuracy = if (ssl_successful[m] > 0) ssl_correct[m] / ssl_successful[m] else NA_real_,
        correct = ssl_correct[m],
        status = if (ssl_successful[m] == ssl_attempted[m] && ssl_attempted[m] > 0) "OK" else if (ssl_successful[m] > 0) "PARTIAL" else "FAIL",
        first_error_message = as.character(ssl_first_error[[m]]),
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(out)
}

run_all <- function(dataset_names) {
  all_results <- list()
  for (ds in dataset_names) {
    cat("\n==============================\n")
    cat("Dataset:", ds, "\n")
    cat("==============================\n")
    res <- tryCatch({
      evaluate_dataset(ds)
    }, error = function(e) {
      data.frame(
        dataset = ds,
        method = NA_character_,
        setting = NA_character_,
        attempted_folds = 0L,
        successful_folds = 0L,
        failed_folds = 0L,
        accuracy = NA_real_,
        correct = 0L,
        status = "DATASET_FAIL",
        first_error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    print(res)
    all_results[[length(all_results) + 1]] <- res
  }
  dplyr::bind_rows(all_results)
}

results_full <- run_all(DATASETS)

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_results <- file.path(OUTPUT_DIR, paste0("full_dd_benchmark_results_", ts, ".csv"))
write.csv(results_full, out_results, row.names = FALSE)
cat("Saved full results to:", out_results, "\n")

summary_methods <- results_full |>
  dplyr::filter(!is.na(method), !is.na(setting)) |>
  dplyr::group_by(setting, method) |>
  dplyr::summarise(
    mean_accuracy_successful = mean(accuracy, na.rm = TRUE),
    successful_datasets = sum(!is.na(accuracy)),
    ok_datasets = sum(status == "OK"),
    partial_datasets = sum(status == "PARTIAL"),
    failed_datasets = sum(status == "FAIL"),
    .groups = "drop"
  ) |>
  dplyr::arrange(setting, failed_datasets, dplyr::desc(mean_accuracy_successful))

print(summary_methods)
out_summary <- file.path(OUTPUT_DIR, paste0("full_dd_benchmark_summary_", ts, ".csv"))
write.csv(summary_methods, out_summary, row.names = FALSE)
cat("Saved method summary to:", out_summary, "\n")

summary_datasets <- results_full |>
  dplyr::filter(!is.na(dataset), !is.na(method)) |>
  dplyr::group_by(dataset) |>
  dplyr::summarise(
    methods_ok = sum(status == "OK"),
    methods_partial = sum(status == "PARTIAL"),
    methods_fail = sum(status == "FAIL"),
    .groups = "drop"
  ) |>
  dplyr::arrange(methods_fail, dplyr::desc(methods_ok))

print(summary_datasets)
out_dataset_summary <- file.path(OUTPUT_DIR, paste0("full_dd_benchmark_dataset_summary_", ts, ".csv"))
write.csv(summary_datasets, out_dataset_summary, row.names = FALSE)
cat("Saved dataset summary to:", out_dataset_summary, "\n")

cat("\nRun completed.\n")

writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")
saveRDS(results_full, "outputs/results_full.rds")




datasets <- unique(results_full$dataset)
dataset_info <- lapply(datasets, function(ds_name) {
  
  ds <- prepare_dataset(ds_name)
  
  n <- nrow(ds$X)
  p <- ncol(ds$X)
  
  class_counts <- table(ds$y)
  class_ratio <- round(min(class_counts) / max(class_counts), 3)
  
  data.frame(
    dataset = ds_name,
    n = n,
    p = p,
    class_0 = class_counts[1],
    class_1 = class_counts[2],
    balance_ratio = class_ratio
  )
})

dataset_table <- bind_rows(dataset_info) |>
  arrange(dataset)

write.csv(dataset_table, "outputs/analysis/tables/dataset_overview.csv", row.names = FALSE)

####

representation_summary <- results_full |>
  filter(status == "OK") |>
  mutate(
    representation = case_when(
      str_detect(method, "\\.orig$") ~ "Original",
      str_detect(method, "\\.poly$") ~ "DD + Poly",
      str_detect(method, "\\.DD$") ~ "DD",
      TRUE ~ "Other"
    )
  ) |>
  group_by(representation) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    n = n()
  ) |>
  arrange(desc(mean_accuracy))

print(representation_summary)

representation_summary <- results_full |>
  filter(status == "OK") |>
  mutate(
    representation = case_when(
      str_detect(method, "\\.orig$") ~ "Original",
      str_detect(method, "\\.poly$") ~ "DD + Poly",
      str_detect(method, "\\.DD$") ~ "DD",
      TRUE ~ "Other"
    )
  ) |>
  group_by(representation) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy = sd(accuracy, na.rm = TRUE),
    n = n()
  )
print(representation_summary)
