# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(microViz))
suppressPackageStartupMessages(library(ggplot2))

# Load Data
ps.dna <- readRDS("phyloseq/ps.dna.rds")

# Remove taxa without phylum-level annotation
colSums(is.na(tax_table(ps.dna)))
table(tax_table(ps.dna)[, "Phylum"], exclude = NULL)

ps.phylum <- subset_taxa(ps.dna, !is.na(Phylum))

# Remove chloroplast and mitochondrial taxa (%in% keeps taxa with NA Order/Family)
ps.taxfilt <- subset_taxa(ps.phylum, !Order %in% c("Chloroplast") &
                            !Family %in% c("Mitochondria"))

# Taxonomy-filtered object, before the prevalence filter: used for alpha diversity,
# since richness depends on the rare taxa that the prevalence filter removes
saveRDS(ps.taxfilt, file = "phyloseq/ps.taxfilt.rds")

# Remove taxa present in less than 5% of samples (== 2 out of 40 samples)
## Check taxa prevalence
prevalence_table <- function(ps) {
  prev <- apply(X = otu_table(ps),
                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
  data.frame(Prevalence = prev,
             TotalAbundance = taxa_sums(ps),
             tax_table(ps))
}

prevalence_by_phylum <- function(prevdf, stage) {
  do.call(rbind, lapply(split(prevdf, prevdf$Phylum), function(df1) {
    data.frame(Stage = stage,
               Phylum = df1$Phylum[1],
               n_taxa = nrow(df1),
               mean_prevalence = mean(df1$Prevalence),
               total_prevalence = sum(df1$Prevalence))
  }))
}

prevdf <- prevalence_table(ps.taxfilt)
prev_summary_BF <- prevalence_by_phylum(prevdf, "Before prevalence filter")
prev_summary_BF

png("results/figures/Prevalence_plot_BF.png", width = 3500, height = 3000, res = 300)
print(ggplot(prevdf, aes(TotalAbundance, Prevalence/nsamples(ps.taxfilt), color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  ggtitle("Prevalence Before Filtering") +
  facet_wrap(~Phylum) + theme(legend.position= "none") +
  theme(axis.title.y = element_text(size = 14, margin = margin(r = 15)),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, margin = margin(b = 20))))
dev.off()

## Filter out taxa prevalence < 5% of samples
ps.taxfilt.prevfilt <- tax_filter(ps.taxfilt, min_prevalence = 0.05)

message("Taxa: ", ntaxa(ps.dna), " -> without NA phylum: ", ntaxa(ps.phylum),
        " -> without chloroplast/mitochondria: ", ntaxa(ps.taxfilt),
        " -> prevalence >= 5%: ", ntaxa(ps.taxfilt.prevfilt))

prevdf_filt <- prevalence_table(ps.taxfilt.prevfilt)
prev_summary_AF <- prevalence_by_phylum(prevdf_filt, "After prevalence filter")
prev_summary_AF

write.csv(rbind(prev_summary_BF, prev_summary_AF),
          "results/tables/Prevalence_by_phylum.csv", quote = F, row.names = F)

png("results/figures/Prevalence_plot_AF.png", width = 3500, height = 3000, res = 300)
print(ggplot(prevdf_filt, aes(TotalAbundance, Prevalence/nsamples(ps.taxfilt.prevfilt), color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  ggtitle("Prevalence After Filtering") +
  facet_wrap(~Phylum) + theme(legend.position= "none") +
  theme(axis.title.y = element_text(size = 14, margin = margin(r = 15)),
        axis.title.x = element_text(size = 14, margin = margin(t = 15)),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        title = element_text(size = 16, margin = margin(b = 20))))
dev.off()

# Save processed object (taxonomy + prevalence filtered): used by the downstream analyses
saveRDS(ps.taxfilt.prevfilt, file = "phyloseq/ps.taxfilt.prevfilt.rds")

# Session information (end of the log)
print(sessionInfo())
