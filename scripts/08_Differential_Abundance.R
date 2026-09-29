# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load package
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(mia))
suppressPackageStartupMessages(library(ANCOMBC))
suppressPackageStartupMessages(library(ggplot2))

# Load files
ps_filt <- readRDS("phyloseq/ps.taxfilt.prevfilt.rds")

# Differential Abundance - ANCOMBC2 (lfc columns are natural-log fold changes, see ?ancombc2)
tse = mia::convertFromPhyloseq(ps_filt)
tse$Group = factor(tse$Group, levels = c("T0", "T12"))

# prv_cut = 0.20: same genus prevalence filter as the beta diversity analysis (present in >= 20% of samples).
# neg_lb = FALSE: recommended by the ANCOM-BC2 documentation for groups with < 30 samples; with TRUE,
# many genera are classified as structural zeros and silently removed from testing.
set.seed(100)
output = ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "Group", rand_formula = "(1 | Patient)",
                  p_adj_method = "BH",
                  prv_cut = 0.20, lib_cut = 1000, s0_perc = 0.05,
                  group = "Group", struc_zero = TRUE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 4, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20,
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

## Structural zeros on genus level: genera absent from ALL samples of one timepoint.
## They cannot be tested for abundance and are reported separately, with their prevalence.
ps_genus <- tax_glom(ps_filt, "Genus", NArm = TRUE)
genus_counts <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) genus_counts <- t(genus_counts)
colnames(genus_counts) <- as.character(tax_table(ps_genus)[, "Genus"])
genus_group <- as.character(sample_data(ps_genus)$Group)

tab_zero = output$zero_ind
write.csv(tab_zero, "results/tables/ANCOMBC2_structural_zeros_all.csv", quote = TRUE, row.names = F)

tab_zero = tab_zero[apply(tab_zero[, -1, drop = FALSE], 1, any), , drop = FALSE]
colnames(tab_zero) <- c("taxon", "absent_in_all_T0", "absent_in_all_T12")
tab_zero$n_present_T0  <- colSums(genus_counts[genus_group == "T0", tab_zero$taxon, drop = FALSE] > 0)
tab_zero$n_present_T12 <- colSums(genus_counts[genus_group == "T12", tab_zero$taxon, drop = FALSE] > 0)
tab_zero
write.csv(tab_zero, "results/tables/ANCOMBC2_structural_zeros.csv", quote = F, row.names = F)

## ANCOM-BC2 primary results (all tested genera)
res_prim = output$res
write.csv(res_prim, "results/tables/ANCOMBC2_results.csv", quote = F, row.names = F)

## Significant taxa
## "Unknown" pools all ASVs without a genus assignment, so it is not interpreted as a genus.
## passed_ss: significance (q < alpha) does not change with the pseudo-count added to zeros
## (0, 0.1, 0.5, 1). For nominal hits (q > alpha at every pseudo-count) it is TRUE; it excludes
## genera whose FDR significance depends on the pseudo-count (artifacts of sparse zeros).
res_genera <- res_prim %>%
  dplyr::filter(taxon != "Unknown")

sig_taxons_padj <- res_genera %>%
  dplyr::filter(`q_GroupT12` < 0.05) %>%
  dplyr::filter(`passed_ss_GroupT12` == TRUE)

## Exploratory: nominal p < 0.05 that also passed the pseudo-count sensitivity analysis
sig_taxons <- res_genera %>%
  dplyr::filter(p_GroupT12 < 0.05) %>%
  dplyr::filter(`passed_ss_GroupT12` == TRUE) %>%
  dplyr::arrange(lfc_GroupT12) %>%
  dplyr::mutate(taxon = factor(taxon, levels = taxon),
                fdr_significant = q_GroupT12 < 0.05)

message("Genera tested: ", nrow(res_genera), "; significant after FDR (robust): ", nrow(sig_taxons_padj),
        "; nominal p < 0.05 and passed_ss: ", nrow(sig_taxons))
write.csv(sig_taxons, "results/tables/ANCOMBC2_nominal_genera.csv", quote = F, row.names = F)

## Plot
png("results/figures/Diff_abundance_ANCOMBC2.png", width = 3000, height = 3000, res = 300)
print(ggplot(sig_taxons, aes(x = lfc_GroupT12, y = taxon,
                       fill = ifelse(lfc_GroupT12 >= 0, "Positive", "Negative"))) +
  geom_bar(stat = "identity", width = 0.7, color = "black",
           position = position_dodge(width = 0.4)) +
  geom_text(data = dplyr::filter(sig_taxons, fdr_significant), aes(label = "*"),
            hjust = ifelse(dplyr::filter(sig_taxons, fdr_significant)$lfc_GroupT12 >= 0, -0.3, 1.3),
            vjust = 0.75, size = 9) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(x = "LFC", y = "Genus",
       title = "T12 relative to T0",
       caption = if (any(sig_taxons$fdr_significant)) "* q < 0.05 (FDR)" else NULL) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Positive" = "purple", "Negative" = "pink"),
                    name = "LFC") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 14, face = "italic"),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 18)))
dev.off()

# Session information (end of the log)
print(sessionInfo())
