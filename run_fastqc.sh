#!/usr/bin/env bash
# Quality control of raw reads: FastQC per sample + MultiQC summary.
# Usage (from anywhere): bash 16SrRNA_analysis/run_fastqc.sh
# Optional: THREADS=8 FASTQC=/path/to/fastqc MULTIQC=/path/to/multiqc
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

raw_dir="raw_fastq"
out_dir="QC_reports"
threads="${THREADS:-12}"
fastqc="${FASTQC:-$(command -v fastqc || echo "$HOME/Downloads/FastQC/fastqc")}"
multiqc="${MULTIQC:-$(command -v multiqc || true)}"

[[ -x "$fastqc" ]] || { echo "FastQC not found; set FASTQC=/path/to/fastqc" >&2; exit 1; }
[[ -n "$multiqc" ]] || { echo "MultiQC not found; set MULTIQC=/path/to/multiqc" >&2; exit 1; }

shopt -s nullglob
fastqs=("$raw_dir"/*.fastq "$raw_dir"/*.fastq.gz)
(( ${#fastqs[@]} )) || { echo "No FASTQ files in $raw_dir" >&2; exit 1; }

mkdir -p "$out_dir"
"$fastqc" --version | tee "$out_dir/fastqc_version.txt"
"$multiqc" --version | tee "$out_dir/multiqc_version.txt"

echo "Running FastQC on ${#fastqs[@]} files with $threads threads"
"$fastqc" --threads "$threads" --outdir "$out_dir" "${fastqs[@]}"

echo "Running MultiQC in $out_dir"
"$multiqc" --force --outdir "$out_dir" "$out_dir"

echo "Done. Report: $project_dir/$out_dir/multiqc_report.html"
