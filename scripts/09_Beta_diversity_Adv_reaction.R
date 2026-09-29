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
suppressPackageStartupMessages(library(zCompositions))
suppressPackageStartupMessages(library(easyCODA))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ecodist))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(permute))

# Load files
ps_filt <- readRDS("phyloseq/ps.taxfilt.prevfilt.rds")

n_perm <- 9999

# Prepare data
## Remove ENL patients (both samples of each patient)
enl_patients <- unique(sample_data(ps_filt)$Patient[sample_data(ps_filt)$Adverse.reaction == "ENL"])
ps <- prune_samples(!sample_data(ps_filt)$Patient %in% enl_patients, ps_filt)
message("ENL patients removed: ", paste(enl_patients, collapse = ", "))

## Agglomerate to genus level
ps_filt_genus <- tax_glom(ps, "Genus", NArm = TRUE)
table <- as.data.frame(otu_table(ps_filt_genus))

## Same zero handling and CLR as the beta diversity analysis (07): keep genera present in >= 20% of
## the samples, replace zeros once using all samples (both timepoints), no sample removed
max_zero_fraction <- 0.8
table <- table[, colMeans(table == 0) <= max_zero_fraction]
table_no0 <- cmultRepl(table, output = "p-counts", z.delete = FALSE, z.warning = 1)
clr_all <- CLR(table_no0, weight = FALSE)$LR  # standard (unweighted) CLR

## Baseline (T0) samples: one row per patient
sample_meta <- as(sample_data(ps_filt_genus), "data.frame")[rownames(clr_all), ]
meta_data <- sample_meta[sample_meta$Group == "T0", ]
clr_pre <- clr_all[meta_data$ID, , drop = FALSE]
rownames(clr_pre) <- meta_data$Patient
rownames(meta_data) <- meta_data$Patient

## Adverse reaction as recorded at T0 (the status at the time of the baseline sample)
meta_data$Adverse.reaction <- factor(meta_data$Adverse.reaction, levels = c("N", "RR"),
                                     labels = c("No reaction", "Reverse reaction"))
meta_data$Run <- factor(meta_data$Run, levels = c("1", "2"))
meta_data$Sex <- factor(meta_data$Sex, levels = c("M", "W"))

message("Patients: ", nrow(meta_data), " (", paste(names(table(meta_data$Adverse.reaction)),
        table(meta_data$Adverse.reaction), sep = " = ", collapse = "; "), "); genera: ", ncol(clr_all))

# Baseline analysis (T0 only)
baseline_dist <- vegdist(clr_pre, method = "euclidean")

## Balance of sequencing run and sex between reaction groups
run_table <- table(Run = meta_data$Run, Reaction = meta_data$Adverse.reaction)
sex_table <- table(Sex = meta_data$Sex, Reaction = meta_data$Adverse.reaction)
fisher_run <- fisher.test(run_table)
fisher_sex <- fisher.test(sex_table)

## Main model: reaction, with permutations restricted within sequencing run (controls the batch
## without spending degrees of freedom; small groups)
set.seed(100)
adonis_reaction_base <- adonis2(baseline_dist ~ Adverse.reaction, data = meta_data,
                                permutations = how(nperm = n_perm, blocks = meta_data$Run))

## Secondary model: reaction adjusted for run and sex (marginal effects)
set.seed(100)
adonis_reaction_base_adj <- adonis2(baseline_dist ~ Adverse.reaction + Run + Sex, data = meta_data,
                                    permutations = n_perm, by = "margin")

set.seed(100)
disp_reaction_base <- permutest(betadisper(baseline_dist, meta_data$Adverse.reaction), permutations = n_perm)
set.seed(100)
disp_batch_base <- permutest(betadisper(baseline_dist, meta_data$Run), permutations = n_perm)
set.seed(100)
disp_sex_base <- permutest(betadisper(baseline_dist, meta_data$Sex), permutations = n_perm)

## Beta diversity
euclidean_pca_baseline <- ecodist::pco(baseline_dist)
axis_labels_baseline <- sprintf("PC%d (%.1f%%)", 1:2, 100 * euclidean_pca_baseline$values[1:2] /
                                  sum(euclidean_pca_baseline$values))

euclidean_pca_baseline_df <- as.data.frame(euclidean_pca_baseline$vectors) %>%
  tibble::rownames_to_column("Patient_ID") %>%
  dplyr::left_join(meta_data, by = c("Patient_ID" = "Patient")) %>%
  dplyr::rename(Patient = Patient_ID)

euclidean_pca_baseline_df$Run <- factor(euclidean_pca_baseline_df$Run, levels = c("1", "2"),
                                        labels = c("Batch 1", "Batch 2"))

## Plot
png("results/figures/Beta_diversity_Adv_reaction_baseline.png", width = 2500, height = 1500, res = 300)
print(ggplot(data = euclidean_pca_baseline_df, aes(x = X1, y = X2)) +
  geom_point(aes(color = Adverse.reaction, shape = Run), size = 3) +
  labs(x = axis_labels_baseline[1], y = axis_labels_baseline[2]) +
  scale_color_manual(name = "Adverse reaction",
                     values = c("No reaction" = "#68228B", "Reverse reaction" = "orange")) +
  scale_shape_discrete(name = "Batch") +
  geom_text_repel(data = euclidean_pca_baseline_df,
                  aes(label = Patient), size = 3.5, max.overlaps = Inf, seed = 100) +
  theme_classic() +
  theme(axis.title = element_text(size = 13, face = "bold", margin = margin(r = 15, t = 15)),
        axis.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11)))
dev.off()

# Save results
sink("results/statistics/Beta_diversity_Adv_reaction_stats.txt")
cat("Patients (ENL excluded: ", paste(enl_patients, collapse = ", "), "): ", nrow(meta_data),
    "; genera (present in >= ", 100 * (1 - max_zero_fraction), "% of samples): ", ncol(clr_all), "\n\n", sep = "")
cat("Sequencing run x reaction:\n")
print(run_table)
cat("Fisher's exact test: p = ", format.pval(fisher_run$p.value, digits = 3), "\n\n", sep = "")
cat("Sex x reaction:\n")
print(sex_table)
cat("Fisher's exact test: p = ", format.pval(fisher_sex$p.value, digits = 3), "\n\n", sep = "")
cat("Baseline Analysis (main): ~ Adverse reaction, permutations within sequencing run\n")
print(adonis_reaction_base)
cat("\n\nBaseline Analysis (secondary): ~ Adverse reaction + Run + Sex, marginal effects\n")
print(adonis_reaction_base_adj)
cat("\n\nHomogeneity Analysis Reaction:\n")
print(disp_reaction_base)
cat("\n\nHomogeneity Analysis Batch:\n")
print(disp_batch_base)
cat("\n\nHomogeneity Analysis Sex:\n")
print(disp_sex_base)
sink()

# Session information (end of the log)
print(sessionInfo())
