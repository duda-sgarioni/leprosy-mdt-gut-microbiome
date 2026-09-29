# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(microbiome))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GGally))

# Load data: taxonomy-filtered object, BEFORE the prevalence filter, because richness
# depends on the rare taxa that the prevalence filter removes
ps_taxfilt <- readRDS("phyloseq/ps.taxfilt.rds")

# Alpha-diversity analysis
# Cassol et al. (2025) offer practical guidelines for selecting alpha diversity indices.
# They classify alpha diversity into four main categories and recommend calculating one metric from each category, using the simplest available metric
# to sufficiently capture within-sample diversity.
# Each index answers a different pre-specified question, so p-values are not adjusted across indices.

## Richness: The number of different microorganisms in a sample
## Dominance (Berger-Parker): relative abundance of the most abundant taxon in the sample (max/total)
## Information: Shannon index is supposed to reflect how many different microbes are in the sample and how evenly they are distributed within a sample
## it is a number that aims to encompass and integrate richness and dominance information
## Phylogenetics: Faith PD is not computed, since no phylogenetic tree was built
## and Faith PD is strongly correlated with the number of observed features
alpha_indices <- function(ps) {
  data.frame(ID = sample_names(ps),
             observed = microbiome::richness(ps, index = "observed")$observed,
             shannon_entropy = microbiome::diversity(ps, index = "shannon")$shannon,
             berger_parker_d = microbiome::dominance(ps, index = "DBP")$dbp)
}

indices <- c("observed", "shannon_entropy", "berger_parker_d")

# Combine different index into one dataframe
alpha_div <- data.frame(as(sample_data(ps_taxfilt), "data.frame"),
                        TotalReads = sample_sums(ps_taxfilt)) %>%
  dplyr::left_join(alpha_indices(ps_taxfilt), by = "ID")

alpha_div$Group <- factor(as.factor(alpha_div$Group), levels = c("T0", "T12"))

# Statistical significance: paired Wilcoxon signed-rank test (T12 vs T0), samples paired by patient.
paired_wilcox <- function(df, index) {
  paired <- df %>%
    dplyr::select(Patient, Group, value = dplyr::all_of(index)) %>%
    tidyr::pivot_wider(names_from = Group, values_from = value)
  stopifnot(!anyNA(paired$T0), !anyNA(paired$T12))
  wilcox.test(paired$T12, paired$T0, paired = TRUE, exact = FALSE)$p.value
}

p_values <- tibble(Index = indices,
                   p_value = sapply(indices, function(i) paired_wilcox(alpha_div, i)))

# Summarised statistics: median (Q1-Q3; min-max) and mean ± SD per timepoint
median_iqr_range <- function(x) {
  q <- quantile(x, c(0.25, 0.75))
  sprintf("%.2f (%.2f-%.2f; %.2f-%.2f)", median(x), q[1], q[2], min(x), max(x))
}
mean_sd <- function(x) sprintf("%.2f ± %.2f", mean(x), sd(x))

alpha_div_summary <- alpha_div %>%
  tidyr::pivot_longer(cols = dplyr::all_of(indices), names_to = "Index", values_to = "value") %>%
  dplyr::group_by(Index, Group) %>%
  dplyr::summarise(n = n(),
                   median_IQR_range = median_iqr_range(value),
                   mean_SD = mean_sd(value),
                   .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(median_IQR_range, mean_SD)) %>%
  dplyr::left_join(p_values, by = "Index") %>%
  dplyr::arrange(match(Index, indices))

write.csv(alpha_div_summary, "results/tables/Alpha_diversity_summary.csv", quote = TRUE, row.names = F)

# Sensitivity analysis: T12 samples are deeper than T0 (see Library_size_metrics.txt), so the indices
# are recomputed after rarefying to the smallest library size (mean of 100 rarefactions)
set.seed(100)
n_rarefactions <- 100
rarefaction_depth <- min(sample_sums(ps_taxfilt))

rarefied_runs <- lapply(seq_len(n_rarefactions), function(i) {
  ps_rare <- rarefy_even_depth(ps_taxfilt, sample.size = rarefaction_depth,
                               rngseed = FALSE, replace = FALSE, trimOTUs = FALSE, verbose = FALSE)
  alpha_indices(ps_rare)
})

alpha_div_rarefied <- alpha_div %>%
  dplyr::select(ID, Patient, Group) %>%
  dplyr::left_join(
    data.frame(ID = rarefied_runs[[1]]$ID,
               sapply(indices, function(i) rowMeans(sapply(rarefied_runs, function(r) r[[i]])))),
    by = "ID")

rarefied_sensitivity <- tibble(Index = indices,
                               p_value_unrarefied = p_values$p_value,
                               p_value_rarefied = sapply(indices, function(i) paired_wilcox(alpha_div_rarefied, i)),
                               rarefaction_depth = rarefaction_depth,
                               n_rarefactions = n_rarefactions)

write.csv(rarefied_sensitivity, "results/tables/Alpha_diversity_rarefied_sensitivity.csv", quote = F, row.names = F)

# diversity_analysis
## Boxplot
index_labels <- c(observed = "Observed Richness",
                  shannon_entropy = "Shannon Entropy",
                  berger_parker_d = "Berger-Parker")

p_labels <- p_values %>%
  dplyr::mutate(label = paste0("p-value = ", signif(p_value, 2)))

plot_alpha <- alpha_div %>%
  tidyr::pivot_longer(cols = dplyr::all_of(indices),
                      names_to = "Index",
                      values_to = "value") %>%
  ggplot(aes(x = Group, y = value, fill = Group)) +
  geom_violin(aes(color = Group, fill = Group), alpha = 0.4, width = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = c("T0" = "#68228B", "T12" = "orange")) +
  scale_color_manual(values = c("T0" = "#68228B", "T12" = "orange")) +
  geom_line(aes(x = Group, y = value, group = Patient), linewidth = 0.7, alpha = 1, color = "grey45") +
  facet_wrap(~ Index, scales = "free_y",
             labeller = as_labeller(index_labels)) +
  labs(x = NULL, y = "Alpha Diversity Measure") +
  geom_text(data = p_labels, aes(label = label), inherit.aes = FALSE,
            x = 1.5, y = Inf, vjust = 1.5, size = 4.3) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        axis.text = element_text(size = 14),
        strip.text = element_text(size=16, face = "bold"),
        panel.grid = element_blank())

png("results/figures/Alpha_diversity.png", width = 4000, height = 3000, res = 300)
print(plot_alpha)
dev.off()

# Correlation between library size and alpha diversity (and among indices)
corr_vars <- c("Library Size" = "TotalReads", "Observed Richness" = "observed",
               "Shannon Entropy" = "shannon_entropy", "Berger-Parker" = "berger_parker_d")

alpha_div_corr <- alpha_div %>%
  dplyr::rename(dplyr::all_of(corr_vars))

corr_pairs <- combn(names(corr_vars), 2, simplify = FALSE)
corr_table <- do.call(rbind, lapply(corr_pairs, function(v) {
  test <- cor.test(alpha_div_corr[[v[1]]], alpha_div_corr[[v[2]]], method = "spearman", exact = FALSE)
  data.frame(Variable_1 = v[1], Variable_2 = v[2], n = nrow(alpha_div_corr),
             spearman_rho = unname(test$estimate), p_value = test$p.value)
}))
write.csv(corr_table, "results/tables/Alpha_diversity_correlations.csv", quote = F, row.names = F)

## Lower panels: points colored by timepoint (visual only), one regression line for all samples
points_by_timepoint <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(aes(color = Group), size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = 0.8) +
    scale_color_manual(name = "Group", values = c("T0" = "#68228B", "T12" = "orange"))
}

## Upper panels: Spearman rho and p-value for all samples (same test as the table above)
spearman_text <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  p_text <- ifelse(test$p.value < 0.001, "< 0.001", sprintf("%.3f", test$p.value))
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 5,
             label = sprintf("\u03c1: %.2f\np-val: %s", test$estimate, p_text)) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
}

corr_plot <- ggpairs(alpha_div_corr, columns = names(corr_vars),
                     upper = list(continuous = spearman_text),
                     diag = list(continuous = "blankDiag"),
                     lower = list(continuous = points_by_timepoint),
                     legend = c(2, 1)) +
  theme(text = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = "bottom")

png("results/figures/Alpha_diversity_correlations.png", width = 4000, height = 3000, res = 300)
print(corr_plot)
dev.off()

# Session information (end of the log)
print(sessionInfo())
