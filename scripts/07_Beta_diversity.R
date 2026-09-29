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

# Agglomerate taxa to genus level
ps_filt_genus <- tax_glom(ps_filt, 'Genus',  NArm = TRUE)

# Beta-diversity - PCA with Aitchison (Euclidean on CLR) distance
table <- as.data.frame(otu_table(ps_filt_genus))

## Keep genera present in >= 20% of the samples (<= 80% zeros): mostly-zero genera would be
## dominated by imputed values. The filter is applied to genera only, so no sample is removed.
max_zero_fraction <- 0.8
table <- table[, colMeans(table == 0) <= max_zero_fraction]
message("Genera kept for CLR: ", ncol(table), " of ", ntaxa(ps_filt_genus))

## Handling zeros (all samples kept: z.delete = FALSE)
table_no0 <- cmultRepl(table, output = "p-counts", z.delete = FALSE, z.warning = 1)

## Extract metadata, in the same sample order as the count table
meta_data <- as(sample_data(ps_filt_genus), "data.frame")[rownames(table_no0), ]
stopifnot(identical(meta_data$ID, rownames(table_no0)))
message("Samples: ", nrow(meta_data), "; patients: ", length(unique(meta_data$Patient)))

## Centered Log Ratio (CLR) transformation (standard, unweighted CLR)
clr_data <- CLR(table_no0, weight = FALSE)

## Compute distances
clr_pca <- clr_data$LR
euclidean_dist <- vegan::vegdist(clr_pca, method = "euclidean")
euclidean_pca <- ecodist::pco(euclidean_dist)

eigenvalues <- euclidean_pca$values
var_explained <- 100 * eigenvalues / sum(eigenvalues)

axis_labels <- sprintf("PC%d (%.1f%%)", seq_along(var_explained),
                       var_explained)
head(round(var_explained, 2))

euclidean_pca_df <- as.data.frame(euclidean_pca$vectors) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(meta_data %>% dplyr::select(ID, Run, Group, Patient), by = "ID")

## Plot
euclidean_pca_df$Group <- factor(as.factor(euclidean_pca_df$Group), levels = c("T0", "T12"))

euclidean_pca_df$Run <- factor(euclidean_pca_df$Run, levels = c("1", "2"),
                                 labels = c("Batch 1", "Batch 2"))

png("results/figures/Beta_diversity.png", width = 3500, height = 2000, res = 300)
print(ggplot(data = euclidean_pca_df,
       aes(x = X1, y = X2)) +
  geom_point(aes(color = Group, shape = Run), size = 3) +
  scale_color_manual(name = "Group",
                     values = c("T0" = "#68228B",
                                "T12" = "orange")) +
  scale_shape_discrete(name = "Batch") +
  geom_text_repel(data = euclidean_pca_df %>% dplyr::filter(Group == "T12"),
                         aes(label = Patient), size = 3.5, max.overlaps = 20, seed = 100) +
  geom_path(aes(group = Patient), linetype = "dashed", alpha = 0.5) +
  labs(x = axis_labels[1], y = axis_labels[2]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, face = "bold", margin = margin(r = 15, t = 15)),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)))
dev.off()

# Statistical significance
meta_data$Group <- factor(as.factor(meta_data$Group), levels = c("T0", "T12"))

meta_data$Run <- factor(meta_data$Run, levels = c("1", "2"),
                        labels = c("Batch 1", "Batch 2"))

meta_data$Sex <- factor(meta_data$Sex, levels = c("M", "W"))

h1 <- how(nperm = n_perm, blocks = meta_data$Patient) # permutation block for paired samples

## Permanova - paired samples (repeated measures: Patient in the model, T0/T12 permuted within patient).
## Only the Group row is interpreted; the Patient p-value is not meaningful under within-patient permutations.
set.seed(100)
adonis_stats <- adonis2(euclidean_dist ~ Patient + Group, data = meta_data, permutations = h1, by = "terms")

set.seed(100)
disp <- permutest(betadisper(euclidean_dist, meta_data$Group), permutations = h1) # check homogeneity

## Partial Redundancy Analysis - paired samples
# Condition(Patient) removes the differences between patients before testing Group, i.e. all effects
# that do not vary within a patient (sex, sequencing run, ...) are separated from the Group effect
pRDA <- rda(clr_pca ~ Group + Condition(Patient), data = meta_data)
pRDA_summ <- summary(pRDA)
pRDA_rsqr <- RsquareAdj(pRDA)$adj.r.squared
set.seed(100)
pRDA_annova <- anova(pRDA, permutations = h1)

## Sensitivity analysis: batch. Run is constant within patient except for the patients whose
## T0 and T12 samples were sequenced in different runs; only there Run can be confounded with Group.
run_per_patient <- tapply(meta_data$Run, meta_data$Patient, function(x) length(unique(x)))
patients_mixed_runs <- names(run_per_patient)[run_per_patient > 1]
keep_same_run <- !meta_data$Patient %in% patients_mixed_runs
meta_same_run <- droplevels(meta_data[keep_same_run, ])
dist_same_run <- as.dist(as.matrix(euclidean_dist)[keep_same_run, keep_same_run])

set.seed(100)
h_same_run <- how(nperm = n_perm, blocks = meta_same_run$Patient)
adonis_same_run <- adonis2(dist_same_run ~ Patient + Group, data = meta_same_run,
                           permutations = h_same_run, by = "terms")

## Sensitivity analysis: replacement treatment (ATS). Does the T0 -> T12 shift appear in patients with
## and without ATS, and does the change itself differ between them? ATS overlaps with adverse reactions
## (reactions can motivate the replacement), so a difference by ATS cannot be attributed to the drugs alone.
### 1. Same paired PERMANOVA within each ATS subgroup (compare effect sizes; power is low in each subgroup)
adonis_by_ats <- lapply(c(Without_ATS = "N", With_ATS = "Y"), function(ats) {
  keep_ats <- meta_data$ATS == ats
  meta_ats <- droplevels(meta_data[keep_ats, ])
  dist_ats <- as.dist(as.matrix(euclidean_dist)[keep_ats, keep_ats])
  set.seed(100)
  adonis2(dist_ats ~ Patient + Group, data = meta_ats,
          permutations = how(nperm = n_perm, blocks = meta_ats$Patient), by = "terms")
})

### 2. Per-patient change (CLR T12 - T0) compared between ATS groups
meta_t0  <- meta_data[meta_data$Group == "T0", ]
meta_t12 <- meta_data[meta_data$Group == "T12", ]
meta_t12 <- meta_t12[match(meta_t0$Patient, meta_t12$Patient), ]
stopifnot(!anyNA(meta_t12$ID))

delta_clr <- clr_pca[meta_t12$ID, , drop = FALSE] - clr_pca[meta_t0$ID, , drop = FALSE]
rownames(delta_clr) <- meta_t0$Patient
delta_dist <- vegan::vegdist(delta_clr, method = "euclidean")
meta_delta <- data.frame(Patient = meta_t0$Patient, ATS = factor(meta_t0$ATS, levels = c("N", "Y")))

set.seed(100)
adonis_ats_delta <- adonis2(delta_dist ~ ATS, data = meta_delta, permutations = n_perm)
set.seed(100)
disp_ats_delta <- permutest(betadisper(delta_dist, meta_delta$ATS), permutations = n_perm)

ats_reaction <- table(ATS = meta_t0$ATS, Adverse_reaction_up_to_T12 = meta_t12$Adverse.reaction)

## Save important info
sink("results/statistics/Beta_diversity_stats.txt")
cat("Samples: ", nrow(meta_data), "; patients: ", length(unique(meta_data$Patient)),
    "; genera (present in >= ", 100 * (1 - max_zero_fraction), "% of samples): ", ncol(table), "\n\n", sep = "")
cat("PERMANOVA (Aitchison distance, ~ Patient + Group, permutations within patient):\n")
print(adonis_stats)
cat("\n\nHomogeneity for Group:\n")
print(disp)
cat("\n\nPartial Redundancy Analysis:\n")
cat("Model summary:\n")
print(pRDA_summ)
cat("Explanatory power:\n")
print(pRDA_rsqr)
cat("\nSignificance test:\n")
print(pRDA_annova)
cat("\n\nSensitivity analysis - excluding patients with T0 and T12 in different sequencing runs (",
    paste(patients_mixed_runs, collapse = ", "), "); patients: ", length(unique(meta_same_run$Patient)), "\n", sep = "")
print(adonis_same_run)
cat("\n\nSensitivity analysis - replacement treatment (ATS)\n")
cat("ATS x adverse reaction (patients):\n")
print(ats_reaction)
for (subgroup in names(adonis_by_ats)) {
  cat("\nPaired PERMANOVA within subgroup ", subgroup, " (patients: ",
      sum(meta_t0$ATS == ifelse(subgroup == "With_ATS", "Y", "N")), "):\n", sep = "")
  print(adonis_by_ats[[subgroup]])
}
cat("\nPERMANOVA of the per-patient change (CLR T12 - T0) ~ ATS:\n")
print(adonis_ats_delta)
cat("\nHomogeneity of the change between ATS groups:\n")
print(disp_ats_delta)
sink()

# Session information (end of the log)
print(sessionInfo())
