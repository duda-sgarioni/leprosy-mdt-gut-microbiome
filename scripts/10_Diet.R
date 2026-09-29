# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Load files
dietary_records <- read.table("metadata/dietary_records.tsv", sep = "\t", header = T, check.names = FALSE)

# Prepare data
## Diet: mean of the record days per patient and timepoint (-1 = T0, -2 = T12).
## MD04, MD11, MD21 and MD24 are excluded: at least one of their timepoints has no dietary record.
mean_diet_long <- dietary_records %>%
  dplyr::select(-matches("^MD04|^MD11|^MD21|^MD24")) %>%
  tidyr::pivot_longer(cols = -c(Record_day, Category),
                      names_to = "sample",
                      values_to = "value") %>%
  dplyr::mutate(id = sub("-[12]$", "", sample),
                time = ifelse(grepl("-1$", sample), "T0", "T12")) %>%
  dplyr::group_by(id, Category, time) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE),
                   .groups = "drop")

mean_diet_long_paired <- mean_diet_long %>%
  tidyr::pivot_wider(names_from = time,
                     values_from = mean_value)

# Group comparison: paired Wilcoxon signed-rank test (T12 vs T0) for every dietary category.
median_iqr_range <- function(x) {
  q <- quantile(x, c(0.25, 0.75))
  sprintf("%.2f (%.2f-%.2f; %.2f-%.2f)", median(x), q[1], q[2], min(x), max(x))
}

categories_list <- unique(mean_diet_long_paired$Category)

summary_diet <- lapply(categories_list, function(categ) {

  data <- mean_diet_long_paired %>% dplyr::filter(Category == categ)
  stopifnot(!anyNA(data$T0), !anyNA(data$T12))
  test <- wilcox.test(data$T12, data$T0, paired = TRUE, exact = FALSE)

  data.frame(Category = categ,
             n = nrow(data),
             T0_median_IQR_range = median_iqr_range(data$T0),
             T12_median_IQR_range = median_iqr_range(data$T12),
             p_value = test$p.value,
             check.names = FALSE)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(p_value = signif(p_value, 3))

print(summary_diet, row.names = FALSE)
write.csv(summary_diet, "results/tables/Diet_summary.csv", quote = TRUE, row.names = FALSE)

# Session information (end of the log)
print(sessionInfo())
