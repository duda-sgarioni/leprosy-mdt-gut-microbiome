# Run from the project root (parent of scripts/) when launched via Rscript.
# When sourcing interactively, set the working directory to the project root first.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg)) {
  setwd(file.path(dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))), ".."))
}
message("Working directory: ", getwd())

# Load Data
meta_data <- read.csv("metadata/metadata_samples.csv")

# One row per patient. Sex, clinical form, age, ATS and OM do not change between timepoints, so they
# are taken from the T0 samples. Adverse reaction is taken from the T12 samples: it describes reactions
# developed up to month 12 of treatment (P17 had no reaction at T0 and developed RR during treatment).
meta_data_T0 <- meta_data[meta_data$Group == "T0", ]
meta_data_T12 <- meta_data[meta_data$Group == "T12", ]
meta_data_T12 <- meta_data_T12[match(meta_data_T0$Patient, meta_data_T12$Patient), ]
stopifnot(!anyNA(meta_data_T12$ID))

n_patients <- nrow(meta_data_T0)

count_percent <- function(x, variable) {
  tab <- table(x)
  data.frame(Variable = variable,
             Category = names(tab),
             Value = sprintf("%d (%.1f%%)", as.integer(tab), 100 * as.numeric(tab) / length(x)))
}

## Age: from the patients' medical records
age <- meta_data_T0$Age
age_shapiro <- shapiro.test(age)
age_quartiles <- quantile(age, c(0.25, 0.75))

clinical_table <- rbind(
  data.frame(Variable = "Patients", Category = "", Value = as.character(n_patients)),
  count_percent(meta_data_T0$Sex, "Sex"),
  data.frame(Variable = "Age (years)", Category = c("mean ± SD", "median (IQR; min-max)"),
             Value = c(sprintf("%.1f ± %.1f", mean(age), sd(age)),
                       sprintf("%.1f (%.1f-%.1f; %d-%d)", median(age), age_quartiles[1], age_quartiles[2],
                               min(age), max(age)))),
  count_percent(meta_data_T0$Leprosy, "Clinical form"),
  count_percent(meta_data_T12$Adverse.reaction, "Adverse reaction (up to T12)"),
  count_percent(meta_data_T0$ATS, "Replacement treatment (ATS)"),
  count_percent(meta_data_T0$OM, "Other medications (OM)"))

print(clinical_table, row.names = FALSE)
cat("\nAge normality (Shapiro-Wilk): W = ", round(age_shapiro$statistic, 3),
    ", p = ", signif(age_shapiro$p.value, 3), "\n", sep = "")

write.csv(clinical_table, "results/tables/Clinical_characteristics.csv", quote = TRUE, row.names = FALSE)

# Session information (end of the log)
print(sessionInfo())
