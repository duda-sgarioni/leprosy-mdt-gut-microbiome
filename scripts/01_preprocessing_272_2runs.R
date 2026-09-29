# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("--check", "--skip-phyloseq", "--phyloseq-only")
if (any(!args %in% allowed_args)) stop("Unknown argument: ", paste(setdiff(args, allowed_args), collapse = ", "))
skip_phyloseq <- "--skip-phyloseq" %in% args
phyloseq_only <- "--phyloseq-only" %in% args
threads <- suppressWarnings(as.integer(Sys.getenv("DADA2_THREADS", "12")))
available_cores <- parallel::detectCores()
if (is.na(threads) || threads < 1L) stop("DADA2_THREADS must be a positive integer")
if (!is.na(available_cores)) threads <- min(threads, available_cores)
options(mc.cores = threads, Ncpus = threads)
message("Using ", threads, " CPU threads; working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(Biostrings))

path <- "raw_fastq/"
filt_path <- "filt_fastq/"
plots_dir <- "dada2_plots/"
phyloseq_dir <- "phyloseq/"
results_dir <- "results/"
figures_dir <- "results/figures"
tables_dir <- "results/tables"
statistics_dir <- "results/statistics"

## Create directories
dirs <- c(filt_path, plots_dir, phyloseq_dir, results_dir, figures_dir, tables_dir, statistics_dir)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

list.files(path)

fnFs <- sort(list.files(path, pattern="\\.fastq(\\.gz)?$", full.names = TRUE))
if (!length(fnFs)) stop("No FASTQ files found in ", path)
fastq_names <- sub("\\.gz$", "", basename(fnFs))

samples_run2 <- c("MD01-1.fastq", "MD01-2.fastq", "MD12-1.fastq", "MD12-2.fastq", "MD14-1.fastq", 
                  "MD17-2.fastq", "MD21-1.fastq", "MD21-2.fastq", "MD22-1.fastq", "MD22-2.fastq",
                  "MD23-1.fastq", "MD23-2.fastq", "MD24-1.fastq", "MD24-2.fastq")

fnFs_run1 <- fnFs[!fastq_names %in% samples_run2]
sample_names_run1 <- sub("\\.fastq(\\.gz)?$", "", basename(fnFs_run1))

fnFs_run2 <- fnFs[fastq_names %in% samples_run2]
sample_names_run2 <- sub("\\.fastq(\\.gz)?$", "", basename(fnFs_run2))

if (!length(fnFs_run1) || !length(fnFs_run2)) stop("Both sequencing runs must contain FASTQ files")
missing_run2 <- setdiff(samples_run2, fastq_names)
if (length(missing_run2)) stop("Missing run 2 files: ", paste(missing_run2, collapse = ", "))
sample_ids <- c(sample_names_run1, sample_names_run2)
if (anyDuplicated(sample_ids)) stop("Duplicate sample IDs in FASTQ files")
references <- c("taxa/silva_nr99_v138.2_toGenus_trainset.fa.gz",
                "taxa/silva_v138.2_assignSpecies.fa.gz")
if (any(!file.exists(references))) stop("Missing SILVA reference: ", paste(references[!file.exists(references)], collapse = ", "))
metadata_file <- Sys.getenv("DADA2_METADATA", "metadata/metadata_samples.csv")
if (!skip_phyloseq) {
  if (!file.exists(metadata_file)) stop("Missing metadata: ", metadata_file,
    ". Supply DADA2_METADATA or use --skip-phyloseq to save preprocessing and taxonomy first.")
  metadata <- read.csv(metadata_file, header = TRUE, colClasses = c(ID = "character"))
  if (!"ID" %in% names(metadata) || anyNA(metadata$ID) || any(!nzchar(metadata$ID)) || anyDuplicated(metadata$ID)) {
    stop("Metadata must contain a nonempty, unique ID column")
  }
  if (length(setdiff(sample_ids, metadata$ID))) stop("Metadata is missing sample IDs: ", paste(setdiff(sample_ids, metadata$ID), collapse = ", "))
  rownames(metadata) <- metadata$ID
}
message("Preflight: ", length(fnFs_run1), " run 1 samples; ", length(fnFs_run2), " run 2 samples; SILVA references present")
if ("--check" %in% args) {
  print(sessionInfo())
  quit(save = "no", status = 0)
}
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))

if (!phyloseq_only) {
# Inspect read quality profiles
png(paste0(plots_dir, "quality_profiles_run1.png"), width = 4000, height = 3000, res = 300)
print(plotQualityProfile(head(fnFs_run1, 26)))
dev.off()

png(paste0(plots_dir, "quality_profiles_run2.png"), width = 4000, height = 3000, res = 300)
print(plotQualityProfile(head(fnFs_run2, 14)))
dev.off()

# Filter and trim
filtFs_run1 <- file.path(filt_path, paste0(sample_names_run1, "_filt.fastq.gz"))
out_run1 <- filterAndTrim(fnFs_run1, filtFs_run1,
                          truncQ = 2,
                          truncLen = 272, # cut before the 806R primer (read pos. 292-19 = 273)
                          trimLeft= 34,
                          maxN = 0, 
                          maxEE = 2,
                          rm.phix=TRUE, compress=TRUE, multithread=threads)
print(out_run1)
sum(out_run1[, "reads.out"])  # Total reads pós-filtragem

filtFs_run2 <- file.path(filt_path, paste0(sample_names_run2, "_filt.fastq.gz"))
out_run2 <- filterAndTrim(fnFs_run2, filtFs_run2,
                          truncQ = 2,
                          truncLen = 272, # cut before the 806R primer (read pos. 292-19 = 273)
                          trimLeft= 34,
                          maxN = 0, 
                          maxEE = 2,
                          rm.phix=TRUE, compress=TRUE, multithread=threads)
print(out_run2)
sum(out_run2[, "reads.out"])  # Total reads after filtering

# Stop before inference if any sample has no usable reads.
if (any(out_run1[, "reads.out"] == 0) || any(out_run2[, "reads.out"] == 0)) {
  stop("Some samples have zero reads after filtering; inspect filtering output before proceeding")
}
saveRDS(list(run1 = out_run1, run2 = out_run2), file.path(results_dir, "filter_counts.rds"))

# Learn the Error Rates
errF_run1 <- learnErrors(filtFs_run1, multithread=threads, nbases = 5e8, verbose = TRUE, MAX_CONSIST = 20)
png(paste0(plots_dir, "error_plot_run1.png"), width=1000, height=600)
print(plotErrors(errF_run1, nominalQ=TRUE))
dev.off()

dada2:::checkConvergence(errF_run1)
png(paste0(plots_dir, "colsums_err_run1.png"), width=1000, height=600)
plot(colSums(errF_run1$trans))
dev.off()

errF_run2 <- learnErrors(filtFs_run2, multithread=threads, nbases = 5e8, verbose = TRUE, MAX_CONSIST = 20)
png(paste0(plots_dir, "error_plot_run2.png"), width=1000, height=600)
print(plotErrors(errF_run2, nominalQ=TRUE))
dev.off()

dada2:::checkConvergence(errF_run2)
png(paste0(plots_dir, "colsums_err_run2.png"), width=1000, height=600)
plot(colSums(errF_run2$trans))
dev.off()

saveRDS(list(run1 = errF_run1, run2 = errF_run2), file.path(results_dir, "error_models.rds"))

# Sample Inference
dadaFs_run1 <- dada(filtFs_run1, err = errF_run1, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread=threads, pool = FALSE)

dadaFs_run2 <- dada(filtFs_run2, err = errF_run2, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, multithread=threads, pool = FALSE)

names(dadaFs_run1) <- sample_names_run1
names(dadaFs_run2) <- sample_names_run2

# Construct sequence table
seqtab_run1 <- makeSequenceTable(dadaFs_run1)
dim(seqtab_run1)
table(nchar(getSequences(seqtab_run1))) # Inspect distribution of sequence lengths

seqtab_run2 <- makeSequenceTable(dadaFs_run2)
dim(seqtab_run2)
table(nchar(getSequences(seqtab_run2))) # Inspect distribution of sequence lengths

# Merge sequence tables from both runs
seqtab <- mergeSequenceTables(seqtab_run1, seqtab_run2)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE, multithread=threads)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## Save infos
sink("info_preprocessing_272.txt")
cat("DADA2 preprocessing summary\n")

cat("Run 1:\n")
cat("  Samples:", nrow(seqtab_run1), "\n")
cat("  Unique ASVs:", ncol(seqtab_run1), "\n\n")

cat("Run 2:\n")
cat("  Samples:", nrow(seqtab_run2), "\n")
cat("  Unique ASVs:", ncol(seqtab_run2), "\n\n")

cat("Merged sequence table:\n")
cat("  Samples:", nrow(seqtab), "\n")
cat("  Unique ASVs:", ncol(seqtab), "\n\n")

cat("After chimera removal:\n")
cat("  Samples:", nrow(seqtab.nochim), "\n")
cat("  Unique ASVs:", ncol(seqtab.nochim), "\n")
cat("  Chimeric ASVs removed:", ncol(seqtab) - ncol(seqtab.nochim), 
    sprintf("(%.1f%%)", 100 * (ncol(seqtab) - ncol(seqtab.nochim)) / ncol(seqtab)), "\n")
cat("  Non-chimeric reads retained:", round(100 * sum(seqtab.nochim) / sum(seqtab), 2), "%\n")

sink()

# Track reads
getN <- function(x) sum(getUniques(x))
track_run1 <- as.data.frame(cbind(out_run1,
                                  sapply(dadaFs_run1, getN)))
track_run2 <- as.data.frame(cbind(out_run2,
                                  sapply(dadaFs_run2, getN)))

colnames(track_run1) <- colnames(track_run2) <- c("input", "filtered", "denoised")

track <- rbind(track_run1, track_run2)
rownames(track) <- c(names(dadaFs_run1), names(dadaFs_run2))
track$nonchim <- rowSums(seqtab.nochim)[rownames(track)]

write.csv(track, "tracking_table_final_272.csv")

saveRDS(seqtab.nochim, file.path(results_dir, "seqtab_nochim.rds"))

# Assign Taxonomy
set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim, "taxa/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread = threads)
taxa <- addSpecies(taxa, "taxa/silva_v138.2_assignSpecies.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

dim(seqtab.nochim)

saveRDS(list(seqtab.nochim = seqtab.nochim, taxa = taxa),
        file.path(results_dir, "preprocessing.rds"))
write.csv(taxa, file.path(tables_dir, "taxonomy.csv"))
} else {
  checkpoint <- readRDS(file.path(results_dir, "preprocessing.rds"))
  seqtab.nochim <- checkpoint$seqtab.nochim
  taxa <- checkpoint$taxa
}
if (skip_phyloseq) {
  message("Preprocessing and taxonomy complete. Phyloseq skipped; supply metadata and run --phyloseq-only later.")
  print(sessionInfo())
  quit(save = "no", status = 0)
}

# Construct the phyloseq object, requiring metadata for every retained sample.
if (length(setdiff(rownames(seqtab.nochim), metadata$ID))) stop("Metadata does not cover all checkpoint samples")
metadata <- metadata[rownames(seqtab.nochim), , drop = FALSE]

rownames(seqtab.nochim) <- sub("_filt\\.fastq\\.gz$", "", rownames(seqtab.nochim))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

saveRDS(ps, file = paste0(phyloseq_dir, "ps.rds")) # sequence name is the sequence itself

## Rename sequences 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps.dna <- merge_phyloseq(ps, dna)
taxa_names(ps.dna) <- paste0("ASV", seq(ntaxa(ps.dna)))
ps.dna

saveRDS(ps.dna, file = paste0(phyloseq_dir, "ps.dna.rds")) # sequence name is ASV and a number

# Session information (end of the log)
print(sessionInfo())
