# Gut microbiota in leprosy patients during treatment (16S rRNA)

Analysis code and results for a longitudinal 16S rRNA gene sequencing study of the gut
microbiota of 20 leprosy patients, sampled before treatment (T0) and at month 12 of
treatment (T12), for 40 samples in total.

This repository contains the scripts, the execution logs and the resulting tables, statistics
and figures. Raw reads, sample metadata and intermediate R objects are not included (see
[Data not included](#data-not-included)).

## Repository structure

```
.
├── run_all_R452.sh              # runs scripts 00-10 in order and writes the logs
├── run_fastqc.sh                # FastQC + MultiQC on the raw reads
├── scripts/                     # analysis scripts (00-10)
├── QC_reports/                  # FastQC reports per sample + MultiQC summary
├── dada2_plots/                 # read quality profiles and DADA2 error models
├── info_preprocessing_272.txt   # DADA2 summary (ASVs per run, chimera removal)
├── tracking_table_final_272.csv # reads per sample at each DADA2 step
└── results/
    ├── figures/
    ├── tables/
    ├── statistics/
    ├── sessionInfo.txt
    └── logs/
        ├── packages_R452.csv    # R packages and versions used
        ├── conda-explicit.txt   # exact Conda environment
        └── 20260928T013414Z/    # one log per script + status.tsv (exit codes)
```

## Pipeline

| Script | Description | Main outputs |
|---|---|---|
| `00_Clinical_Sociodemographic.R` | Clinical and sociodemographic characteristics of the patients | `tables/Clinical_characteristics.csv` |
| `01_preprocessing_272_2runs.R` | DADA2 per sequencing run (forward reads, `trimLeft = 34`, `truncLen = 272`, `maxEE = 2`), merge of runs, chimera removal, SILVA v138.2 taxonomy, phyloseq objects | `dada2_plots/`, `tables/taxonomy.csv`, `info_preprocessing_272.txt`, `tracking_table_final_272.csv` |
| `02_Prepare_phyloseq.R` | Removal of taxa without phylum assignment, chloroplasts and mitochondria; prevalence filter (≥ 5% of samples) | `figures/Prevalence_plot_*.png`, `tables/Prevalence_by_phylum.csv` |
| `03_Explore_data.R` | Library sizes and rarefaction curves | `statistics/Library_size_metrics.txt`, `figures/Library_size.png`, `figures/Rarefaction_curve.png` |
| `04_Microbiome_overview.R` | Dominant phyla, genera and ASVs and their relative abundance | `tables/Abundance_*.csv`, `tables/Average_abundance_*.csv` |
| `05_Community_composition.R` | Relative abundance barplot at phylum level | `figures/Community_composition_barplot_RA.png` |
| `06_Alpha_diversity.R` | Observed richness, Shannon entropy and Berger-Parker dominance; paired Wilcoxon tests (T12 vs T0), rarefied sensitivity analysis and correlations between indices | `tables/Alpha_diversity_*.csv`, `figures/Alpha_diversity*.png` |
| `07_Beta_diversity.R` | Aitchison distance (CLR, genus level); paired PERMANOVA and dispersion tests, partial RDA, sensitivity analyses for sequencing run and replacement treatment | `statistics/Beta_diversity_stats.txt`, `figures/Beta_diversity.png` |
| `08_Differential_Abundance.R` | ANCOM-BC2 at genus level (T12 vs T0, patient as random effect, BH correction) | `tables/ANCOMBC2_*.csv`, `figures/Diff_abundance_ANCOMBC2.png` |
| `09_Beta_diversity_Adv_reaction.R` | Baseline (T0) beta diversity by adverse reaction during treatment | `statistics/Beta_diversity_Adv_reaction_stats.txt`, `figures/Beta_diversity_Adv_reaction_baseline.png` |
| `10_Diet.R` | Dietary intake at T0 and T12; paired Wilcoxon tests | `tables/Diet_summary.csv` |

Paths in the Outputs column are relative to `results/` unless they are at the repository root.

## Running the analysis

The scripts expect the project root as working directory and the input data described in
[Data not included](#data-not-included). From the folder that contains `16SrRNA_analysis/`, run:

```bash
bash 16SrRNA_analysis/run_all_R452.sh
```

The runner uses the R installation in `../.tools/r-4.5.2` (relative to the project folder), checks
the R version, creates the output folders and runs scripts 00 to 10 in separate processes. It stops
if a step fails. Existing result files may be overwritten; logs go to a new folder for each run.

Logs are written to `results/logs/<UTC date>/`, with one file per script and a `status.tsv` with
start/end times and exit codes (0 means success). DADA2 uses 12 threads by default, which can be
changed with the `DADA2_THREADS` environment variable.

Quality control of the raw reads is run separately with `run_fastqc.sh`, which writes the FastQC
and MultiQC reports to `QC_reports/` (`multiqc_report.html` summarizes all samples).

## Environment

- R 4.5.2.
- Binary packages from conda-forge and bioconda.
- `easyCODA` 0.40.2 and `ellipse` 0.5.0: installed from CRAN.
- `microViz` 0.13.1: installed from <https://david-barnett.r-universe.dev>.
- `rbiom` 2.2.1: installed from the CRAN archive for compatibility with `mia` 1.18.0, which imports
  the `unifrac` function removed in `rbiom` 3.1.0.

`results/logs/packages_R452.csv` lists the packages actually installed.
`results/logs/conda-explicit.txt` lists the packages provided by Conda. To recreate the environment,
the source packages listed above must also be installed, including replacing the Conda `rbiom` with
version 2.2.1.

Some binaries were built with R 4.5.3 and emit a warning when loaded; the scripts run on R 4.5.2.

## Data not included

- **Raw reads** (`raw_fastq/`): to be deposited in a public repository; the accession number will be added here.
- **Sample metadata** (`metadata/`): clinical, genetic and dietary data of the participants are not
  shared publicly.
- **Reference databases** (`taxa/`): SILVA v138.2 formatted for DADA2
  (`silva_nr99_v138.2_toGenus_trainset.fa.gz` and `silva_v138.2_assignSpecies.fa.gz`), available at
  <https://doi.org/10.5281/zenodo.14169026>.
- **Intermediate files**: filtered reads and R objects (`.rds`), which are regenerated by the scripts.
