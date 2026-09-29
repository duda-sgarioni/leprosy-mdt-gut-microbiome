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

# Load files
ps_filt <- readRDS("phyloseq/ps.taxfilt.prevfilt.rds")

# Timepoints: T0 (before treatment), T12 (month 12 of treatment)
timepoint <- as.character(sample_data(ps_filt)$Group)
stopifnot(all(timepoint %in% c("T0", "T12")))

# Dominant taxa and their overall relative abundance
# Relative abundance is computed per sample BEFORE aggregating, so every level is a
# percentage of all reads; reads without an assignment at a rank go to "Unassigned".
# Species level is not reported: exact-match species assignment covers ~10% of reads.
ps.ra <- transform_sample_counts(ps_filt, function(x) x/sum(x))
ra <- as(otu_table(ps.ra), "matrix")
if (taxa_are_rows(ps.ra)) ra <- t(ra)  # samples x taxa
tax <- as.data.frame(as(tax_table(ps.ra), "matrix"))
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

average_abundance <- function(rank) {
  if (rank == "ASV") {
    lineage <- cbind(ASV = colnames(ra), tax[, ranks])
    key <- colnames(ra)
  } else {
    lineage <- tax[, ranks[seq_len(match(rank, ranks))], drop = FALSE]
    key <- ifelse(is.na(tax[[rank]]), "Unassigned", do.call(paste, c(lineage, sep = ";")))
  }
  agg <- t(rowsum(t(ra), key))  # samples x taxa of this rank
  out <- data.frame(lineage[match(colnames(agg), key), , drop = FALSE],
                    Abundance_all = 100 * colMeans(agg),
                    Abundance_T0  = 100 * colMeans(agg[timepoint == "T0", , drop = FALSE]),
                    Abundance_T12 = 100 * colMeans(agg[timepoint == "T12", , drop = FALSE]),
                    row.names = NULL)
  unassigned <- colnames(agg) == "Unassigned"
  out[unassigned, names(lineage)] <- NA
  out[unassigned, rank] <- "Unassigned"
  stopifnot(abs(sum(out$Abundance_all) - 100) < 1e-6)
  # Most abundant first, "Unassigned" last
  out[order(unassigned, -out$Abundance_all), ]
}

df_asv    <- average_abundance("ASV")
df_phylum <- average_abundance("Phylum")
df_genus  <- average_abundance("Genus")

## Summary of the number of taxa above each mean abundance threshold ("Unassigned" not counted)
count_taxa <- function(df, rank) {
  x <- df$Abundance_all[df[[rank]] != "Unassigned"]
  c(length(x), sum(x > 0.01), sum(x > 0.1), sum(x > 1))
}

df.count <- data.frame(Included = c("All", "> 0.01%", "> 0.1%", "> 1%"),
                       Phylum = count_taxa(df_phylum, "Phylum"),
                       Genus = count_taxa(df_genus, "Genus"),
                       ASV = count_taxa(df_asv, "ASV"))
print(df.count)
write.csv(df.count, "results/tables/Abundance_phyla_genera_ASV.csv", quote = F, row.names = F)

# Average abundance (% of all reads): overall, T0 and T12
## Phyla
print(df_phylum %>% dplyr::select(Phylum, starts_with("Abundance")), digits = 3, row.names = FALSE)
write.csv(df_phylum, "results/tables/Average_abundance_phylum.csv", quote = F, row.names = F)

## Genus: top 10 dominating
print(df_genus %>% dplyr::select(Family, Genus, starts_with("Abundance")) %>% dplyr::slice_head(n = 10),
      digits = 3, row.names = FALSE)
cat("Unassigned at genus level (% of reads):",
    round(df_genus$Abundance_all[df_genus$Genus == "Unassigned"], 2), "\n")
write.csv(df_genus, "results/tables/Average_abundance_genus.csv", quote = F, row.names = F)

## ASV
write.csv(df_asv, "results/tables/Average_abundance_ASV.csv", quote = F, row.names = F)

# Session information (end of the log)
print(sessionInfo())
