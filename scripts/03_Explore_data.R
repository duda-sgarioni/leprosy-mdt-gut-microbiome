# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

# Load files
ps_taxfilt <- readRDS("phyloseq/ps.taxfilt.rds")

# Library size
sample_counts = data.table(as(sample_data(ps_taxfilt), "data.frame"),
                           TotalReads = sample_sums(ps_taxfilt), keep.rownames = TRUE)

sample_counts$Group <- factor(as.factor(sample_counts$Group), levels = c("T0", "T12"))

## Summary statistics
depth_summary <- function(x) {
  q <- quantile(x, c(0.25, 0.75))
  sprintf("median (IQR; min-max) = %.0f (%.0f-%.0f; %.0f-%.0f)", median(x), q[1], q[2], min(x), max(x))
}

group_levels <- levels(sample_counts$Group)
depth_t1 <- sample_counts[sample_counts$Group == group_levels[1], ]
depth_t2 <- sample_counts[sample_counts$Group == group_levels[2], ]
depth_t2 <- depth_t2[match(depth_t1$Patient, depth_t2$Patient), ] # pair samples by patient
stopifnot(identical(as.character(depth_t1$Patient), as.character(depth_t2$Patient)))
depth_paired_test <- wilcox.test(depth_t2$TotalReads, depth_t1$TotalReads, paired = TRUE)
depth_run_test <- wilcox.test(TotalReads ~ Run, data = sample_counts)

sink("results/statistics/Library_size_metrics.txt")
cat("Library size (reads per sample)\n\n")
cat("All samples (n = ", nrow(sample_counts), "): ", depth_summary(sample_counts$TotalReads), "\n\n", sep = "")

cat("By group:\n")
for (g in group_levels) {
  x <- sample_counts$TotalReads[sample_counts$Group == g]
  cat("  ", g, " (n = ", length(x), "): ", depth_summary(x), "\n", sep = "")
}
cat("\nPaired Wilcoxon signed-rank test (", group_levels[2], " vs ", group_levels[1], "):\n", sep = "")
cat("  V = ", depth_paired_test$statistic, ", p = ", format.pval(depth_paired_test$p.value, digits = 3), "\n", sep = "")
cat("  ", group_levels[2], " deeper than ", group_levels[1], " in ",
    sum(depth_t2$TotalReads > depth_t1$TotalReads), " of ", nrow(depth_t1), " patients\n\n", sep = "")

cat("By sequencing run:\n")
for (r in sort(unique(sample_counts$Run))) {
  x <- sample_counts$TotalReads[sample_counts$Run == r]
  cat("  Run ", r, " (n = ", length(x), "): ", depth_summary(x), "\n", sep = "")
}
cat("\nWilcoxon rank-sum test (Run 1 vs Run 2): W = ", depth_run_test$statistic,
    ", p = ", format.pval(depth_run_test$p.value, digits = 3), "\n\n", sep = "")
cat("Samples per run and group:\n")
print(table(Run = sample_counts$Run, Group = sample_counts$Group))
sink()

## Plot
png("results/figures/Library_size.png", width = 3000, height = 3000, res = 300)
print(ggplot(sample_counts, aes(x = Group, y = TotalReads, colour = Patient)) + 
  geom_violin(fill = "grey80", color = "grey40") + labs(x = NULL, y = "Library Size") + 
  geom_point(size = 2.5) + 
  geom_line(aes(group = Patient), linewidth = 0.7, alpha = 1) + 
  geom_text_repel(data = sample_counts %>% dplyr::filter(Group == "T12"), 
                  aes(label = Patient), size = 3.5, max.overlaps = 20, seed = 100) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)), 
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)))
dev.off()

# Rarefaction curve (before the prevalence filter, consistent with alpha diversity)
otu <- as(otu_table(ps_taxfilt), "matrix")
raremax <- min(rowSums(otu))

png("results/figures/Rarefaction_curve.png", width = 4000, height = 3000, res = 300)
rarecurve(otu, step = 20, col = "red", cex = 0.7, cex.axis = 1.4, cex.lab = 1.5, label = TRUE)
abline(v = raremax, col = "black", lty = 3, lwd = 2)
dev.off()

# Session information (end of the log)
print(sessionInfo())
