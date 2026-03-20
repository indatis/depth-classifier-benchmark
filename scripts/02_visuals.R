# ============================================================
# 03_analysis_visuals_dd_benchmark.R
# Visual and summary analysis for the DD benchmark
# Uses:
#   - all results (including failures)
#   - successful-only results
# ============================================================

rm(list = ls())

# -------------------------------
# Packages
# -------------------------------
pkg <- c("dplyr", "ggplot2", "tidyr", "readr", "forcats", "reshape2", "pheatmap")
new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
invisible(lapply(pkg, require, character.only = TRUE))

# -------------------------------
# Paths
# -------------------------------
OUTPUT_DIR <- "outputs"
ANALYSIS_DIR <- file.path(OUTPUT_DIR, "analysis")
FIG_DIR <- file.path(ANALYSIS_DIR, "figures")
TAB_DIR <- file.path(ANALYSIS_DIR, "tables")

dirs_to_create <- c(ANALYSIS_DIR, FIG_DIR, TAB_DIR)
for (d in dirs_to_create) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------------
# Load benchmark results
# -------------------------------
# Adjust this if needed
results_file <- file.path(OUTPUT_DIR, "full_dd_benchmark_results_20260320_113144.csv")

results_full <- read.csv(results_file, stringsAsFactors = FALSE)

# -------------------------------
# Basic cleaning
# -------------------------------
results_full <- results_full |>
  mutate(
    method = as.character(method),
    setting = as.character(setting),
    dataset = as.character(dataset),
    status = as.character(status)
  )

# Successful / usable subset for fair accuracy comparison
results_success <- results_full |>
  filter(!is.na(method), !is.na(setting), !is.na(dataset), !is.na(accuracy), status != "FAIL")

# Keep only rows with actual methods
results_methods <- results_full |>
  filter(!is.na(method), !is.na(setting), !is.na(dataset))

# -------------------------------
# Method grouping
# -------------------------------
add_method_group <- function(df) {
  df |>
    mutate(
      method_group = case_when(
        method == "DD.alpha" ~ "DD classifier",
        grepl("orig$", method) ~ "Original representation",
        grepl("\\.DD$", method) ~ "DD features",
        grepl("poly$", method) ~ "DD + polynomial features",
        TRUE ~ "Other"
      )
    )
}

results_full2 <- add_method_group(results_methods)
results_success2 <- add_method_group(results_success)

# -------------------------------
# 1. METHOD SUMMARY TABLES
# -------------------------------

# All results: includes failures
summary_methods_all <- results_full2 |>
  group_by(setting, method) |>
  summarise(
    mean_accuracy_successful = mean(accuracy, na.rm = TRUE),
    median_accuracy_successful = median(accuracy, na.rm = TRUE),
    sd_accuracy_successful = sd(accuracy, na.rm = TRUE),
    successful_datasets = sum(!is.na(accuracy)),
    ok_datasets = sum(status == "OK"),
    partial_datasets = sum(status == "PARTIAL"),
    failed_datasets = sum(status == "FAIL"),
    .groups = "drop"
  ) |>
  arrange(setting, failed_datasets, desc(mean_accuracy_successful))

write.csv(summary_methods_all,
          file.path(TAB_DIR, "summary_methods_all.csv"),
          row.names = FALSE)

# Successful-only results: fair accuracy comparison
summary_methods_success <- results_success2 |>
  group_by(setting, method) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    median_accuracy = median(accuracy, na.rm = TRUE),
    sd_accuracy = sd(accuracy, na.rm = TRUE),
    n_datasets = n_distinct(dataset),
    .groups = "drop"
  ) |>
  arrange(setting, desc(mean_accuracy))

write.csv(summary_methods_success,
          file.path(TAB_DIR, "summary_methods_success_only.csv"),
          row.names = FALSE)

# -------------------------------
# 2. DATASET SUMMARY TABLES
# -------------------------------

summary_datasets_all <- results_full2 |>
  group_by(dataset) |>
  summarise(
    methods_ok = sum(status == "OK"),
    methods_partial = sum(status == "PARTIAL"),
    methods_fail = sum(status == "FAIL"),
    .groups = "drop"
  ) |>
  arrange(methods_fail, desc(methods_ok))

write.csv(summary_datasets_all,
          file.path(TAB_DIR, "summary_datasets_all.csv"),
          row.names = FALSE)

summary_datasets_success <- results_success2 |>
  group_by(dataset) |>
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    best_accuracy = max(accuracy, na.rm = TRUE),
    n_successful_methods = n(),
    .groups = "drop"
  ) |>
  arrange(desc(best_accuracy), desc(mean_accuracy))

write.csv(summary_datasets_success,
          file.path(TAB_DIR, "summary_datasets_success_only.csv"),
          row.names = FALSE)

# -------------------------------
# 3. FAILURE TABLE
# -------------------------------
failure_table <- results_full2 |>
  count(method, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0)

write.csv(failure_table,
          file.path(TAB_DIR, "failure_table_by_method.csv"),
          row.names = FALSE)

# -------------------------------
# 4. BEST METHOD PER DATASET
# -------------------------------

# Successful-only best method(s)
best_method_per_dataset <- results_success2 |>
  group_by(dataset) |>
  filter(accuracy == max(accuracy, na.rm = TRUE)) |>
  ungroup()

best_method_counts <- best_method_per_dataset |>
  count(method, sort = TRUE)

write.csv(best_method_per_dataset,
          file.path(TAB_DIR, "best_method_per_dataset.csv"),
          row.names = FALSE)

write.csv(best_method_counts,
          file.path(TAB_DIR, "best_method_counts.csv"),
          row.names = FALSE)

# -------------------------------
# 5. BOXPLOT: ALL SUCCESSFUL RESULTS
# -------------------------------
p_box_success <- ggplot(results_success2,
                        aes(x = forcats::fct_reorder(method, accuracy, .fun = median, .na_rm = TRUE),
                            y = accuracy,
                            fill = setting)) +
  geom_boxplot(outlier.alpha = 0.5) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Accuracy distribution across datasets (successful results only)",
    x = "Method",
    y = "Accuracy"
  )

ggsave(file.path(FIG_DIR, "boxplot_accuracy_success_only.png"),
       p_box_success, width = 10, height = 6, dpi = 300)

# -------------------------------
# 6. BOXPLOT: SUCCESSFUL RESULTS BY METHOD GROUP
# -------------------------------
p_box_group <- ggplot(results_success2,
                      aes(x = method_group, y = accuracy, fill = setting)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Accuracy by representation/classification group",
    x = "Method group",
    y = "Accuracy"
  )

ggsave(file.path(FIG_DIR, "boxplot_accuracy_by_group_success_only.png"),
       p_box_group, width = 9, height = 6, dpi = 300)

# -------------------------------
# 7. BOXPLOT: SUPERVISED VS SEMI-SUPERVISED
# -------------------------------
p_setting <- ggplot(results_success2,
                    aes(x = setting, y = accuracy, fill = setting)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Accuracy by learning setting (successful results only)",
    x = "Setting",
    y = "Accuracy"
  ) +
  guides(fill = "none")

ggsave(file.path(FIG_DIR, "boxplot_setting_success_only.png"),
       p_setting, width = 7, height = 5, dpi = 300)

# -------------------------------
# 8. FAILURE BARPLOT
# -------------------------------
p_fail <- ggplot(results_full2,
                 aes(x = forcats::fct_infreq(method), fill = status)) +
  geom_bar(position = "stack") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Method robustness: OK / PARTIAL / FAIL counts",
    x = "Method",
    y = "Number of datasets"
  )

ggsave(file.path(FIG_DIR, "barplot_method_status_counts.png"),
       p_fail, width = 9, height = 6, dpi = 300)

# -------------------------------
# 9. BEST-METHOD COUNT BARPLOT
# -------------------------------
p_best <- ggplot(best_method_counts,
                 aes(x = forcats::fct_reorder(method, n), y = n)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Number of datasets won by each method",
    x = "Method",
    y = "Win count"
  )

ggsave(file.path(FIG_DIR, "barplot_best_method_counts.png"),
       p_best, width = 8, height = 5, dpi = 300)


# -------------------------------
# 10. HEATMAP: SUCCESSFUL RESULTS ONLY (IMPROVED)
# -------------------------------
library(pheatmap)

heatmap_success_df <- results_success2 |>
  select(dataset, method, accuracy) |>
  pivot_wider(names_from = method, values_from = accuracy)

heatmap_success_mat <- as.matrix(heatmap_success_df[, -1])
rownames(heatmap_success_mat) <- heatmap_success_df$dataset

png(file.path(FIG_DIR, "heatmap_accuracy_success_only.png"),
    width = 1800, height = 2400, res = 220)

pheatmap(
  heatmap_success_mat,
  cluster_rows = T,
  cluster_cols = T,
  color = colorRampPalette(c("darkred", "yellow", "darkgreen"))(100),
  main = "LOOCV accuracy heatmap across datasets and methods",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  na_col = "grey90",
  angle_col = 45
)

dev.off()

# -------------------------------
# 11. HEATMAP: STATUS MATRIX (IMPROVED)
# -------------------------------
status_numeric <- results_full2 |>
  mutate(
    status_code = case_when(
      status == "OK" ~ 2,
      status == "PARTIAL" ~ 1,
      status == "FAIL" ~ 0,
      TRUE ~ NA_real_
    )
  ) |>
  select(dataset, method, status_code) |>
  pivot_wider(names_from = method, values_from = status_code)

status_mat <- as.matrix(status_numeric[, -1])
rownames(status_mat) <- status_numeric$dataset

status_colors <- c("darkred", "orange", "darkgreen")

png(file.path(FIG_DIR, "heatmap_status_matrix.png"),
    width = 1800, height = 2400, res = 220)

pheatmap(
  status_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = status_colors,
  breaks = c(-0.5, 0.5, 1.5, 2.5),
  main = "Execution status heatmap across datasets and methods",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 11,
  border_color = NA,
  legend_breaks = c(0, 1, 2),
  legend_labels = c("FAIL", "PARTIAL", "OK"),
  angle_col = 45,
  na_col = "grey90"
)

dev.off()

# -------------------------------
# 12. METHOD MEAN ACCURACY BARPLOT (SUCCESSFUL ONLY)
# -------------------------------
p_mean_acc <- summary_methods_success |>
  ggplot(aes(x = forcats::fct_reorder(method, mean_accuracy), y = mean_accuracy, fill = setting)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean accuracy by method (successful results only)",
    x = "Method",
    y = "Mean accuracy"
  )

ggsave(file.path(FIG_DIR, "barplot_mean_accuracy_success_only.png"),
       p_mean_acc, width = 9, height = 6, dpi = 300)

# -------------------------------
# 13. OPTIONAL: FRIEDMAN TEST
# -------------------------------
# Uses successful-only results.
# Because some methods failed on some datasets, restrict to methods that succeeded on all datasets
methods_complete <- summary_methods_all |>
  filter(successful_datasets == max(successful_datasets, na.rm = TRUE)) |>
  pull(method)

friedman_df <- results_success2 |>
  filter(method %in% methods_complete) |>
  select(dataset, method, accuracy)

friedman_result <- tryCatch({
  friedman.test(accuracy ~ method | dataset, data = friedman_df)
}, error = function(e) e)

sink(file.path(TAB_DIR, "friedman_test_output.txt"))
print(friedman_result)
sink()

# -------------------------------
# 14. OPTIONAL PCA EXAMPLE FOR ONE DATASET
# -------------------------------
# Change dataset if you want
prepare_dataset <- function(dataset_name) {
  data <- getdata(dataset_name)
  colnames(data)[ncol(data)] <- "C"
  
  if (dataset_name %in% c("cloud", "gemsen_MvsF", "wine_2vs3")) {
    data$C <- as.numeric(data$C)
  }
  
  data <- na.omit(data)
  X <- subset(data, select = -C)
  X <- as.data.frame(lapply(X, as.numeric))
  y <- as.factor(data$C)
  levels(y) <- c("0", "1")
  list(X = X, y = y)
}

ds <- prepare_dataset("bupa")
pca <- prcomp(ds$X, scale. = TRUE)

png(file.path(FIG_DIR, "pca_bupa.png"),
    width = 900, height = 700, res = 140)
plot(
  pca$x[, 1:2],
  col = as.numeric(ds$y),
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA projection of bupa"
)
legend("topright", legend = levels(ds$y), col = 1:length(levels(ds$y)), pch = 19)
dev.off()

# -------------------------------
# 15. SAVE CLEAN RESULTS TOO
# -------------------------------
write.csv(results_full2,
          file.path(TAB_DIR, "results_full_with_groups.csv"),
          row.names = FALSE)

write.csv(results_success2,
          file.path(TAB_DIR, "results_success_only_with_groups.csv"),
          row.names = FALSE)

cat("\nAnalysis completed.\n")
cat("Figures saved in:", FIG_DIR, "\n")
cat("Tables saved in:", TAB_DIR, "\n")


##### Extra
png(file.path(FIG_DIR, "best_method.png"),
    width = 900, height = 700, res = 140)
  ggplot(best_method_counts, aes(x = reorder(method, n), y = n)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    ggtitle("Number of datasets where method is best")
dev.off()

library(dplyr)
compare_df <- results_success |>
  select(dataset, method, accuracy) |>
  tidyr::pivot_wider(names_from = method, values_from = accuracy)
mean(compare_df$DD.alpha > compare_df$LDA.orig, na.rm = TRUE)

#
compare_df$diff_DD_vs_LDA <- compare_df$DD.alpha - compare_df$LDA.orig

png(file.path(FIG_DIR, "acc_difference.png"),
    width = 900, height = 700, res = 140)
ggplot(compare_df, aes(x = diff_DD_vs_LDA)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  ggtitle("DD.alpha vs LDA.orig (accuracy difference)")
dev.off()

png(file.path(FIG_DIR, "acc_difference_scatter.png"),
    width = 900, height = 700, res = 140)
ggplot(compare_df, aes(x = LDA.orig, y = DD.alpha)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_minimal() +
  ggtitle("DD.alpha vs LDA.orig per dataset")
dev.off()

mean(compare_df$DD.alpha - compare_df$LDA.orig, na.rm = TRUE)
median(compare_df$DD.alpha - compare_df$LDA.orig, na.rm = TRUE)


