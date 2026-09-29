# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load packages
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggnewscale))

# Load files
ps_filt <- readRDS("phyloseq/ps.taxfilt.prevfilt.rds")

# Agglomerate taxa to phylum level
ps_phylum <- tax_glom(ps_filt, "Phylum", NArm = TRUE)

# Composition barplot
## Compute relative abundance
ps_phylum_RA <- transform_sample_counts(ps_phylum, function(x) x/sum(x))

df_barplot <- psmelt(ps_phylum_RA)

df_barplot$Group <- factor(as.factor(df_barplot$Group), levels = c("T0", "T12"))

phylum_colors <- c("Actinomycetota" = "#E41A1C",
                   "Bacillota"= "#377EB8",
                   "Bacteroidota" = "#4DAF4A",
                   "Cyanobacteriota" = "#984EA3",
                   "Fusobacteriota" = "#FF7F00",
                   "Methanobacteriota" = "#A65628",
                   "Pseudomonadota" = "#F781BF",
                   "Synergistota" = "#999999",
                   "Thermodesulfobacteriota" = "#66C2A5",
                   "Thermoplasmatota" = "#FFD92F",
                   "Verrucomicrobiota" = "#A6CEE3")

## Every phylum in the data needs a color
missing_colors <- setdiff(unique(df_barplot$Phylum), names(phylum_colors))
if (length(missing_colors)) {
  stop("Add a color to phylum_colors for: ", paste(missing_colors, collapse = ", "))
}

## Plot
png("results/figures/Community_composition_barplot_RA.png", width = 6000, height = 3000, res = 300)
print(ggplot(df_barplot, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = NULL, y = "Relative Abundance") +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  scale_fill_manual(values = phylum_colors,
                    guide = guide_legend(label.theme = element_text(face = "italic", size = 14))) +
  new_scale_fill() +
  geom_tile(data = df_barplot, aes(x = Sample, y = 1.1, fill = Group), height = 0.1) +
  scale_fill_manual(values = c("T0" = "#68228B", "T12" = "orange")) +
  theme_minimal() +
  theme(plot.margin = margin(t = 23, r = 5, b = 5, l = 5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, color = "black", angle = 90, margin = margin(t = 10)),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold")))
dev.off()

# Session information (end of the log)
print(sessionInfo())
