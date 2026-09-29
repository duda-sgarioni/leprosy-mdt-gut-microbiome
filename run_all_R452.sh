#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
runtime_dir="$(dirname -- "$project_dir")/.tools"
export PATH="$runtime_dir/r-4.5.2/bin:$PATH"
export R_LIBS_USER="$runtime_dir/r-4.5.2/lib/R/library"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export DADA2_THREADS="${DADA2_THREADS:-12}"
cd "$project_dir"
Rscript --vanilla -e 'stopifnot(getRversion() == "4.5.2")'
run_id="$(date -u +%Y%m%dT%H%M%SZ)"
log_dir="$project_dir/results/logs/$run_id"
mkdir -p "$log_dir" results/figures results/tables results/statistics
printf 'script\tstart_utc\tend_utc\texit_code\n' > "$log_dir/status.tsv"
for script in scripts/[0-9][0-9]_*.R; do
  name="$(basename "$script" .R)"
  started="$(date -u +%FT%TZ)"
  printf '[%s] Starting %s\n' "$started" "$name"
  code=0
  Rscript --vanilla "$script" > "$log_dir/$name.log" 2>&1 || code=$?
  printf '%s\t%s\t%s\t%s\n' "$name" "$started" "$(date -u +%FT%TZ)" "$code" >> "$log_dir/status.tsv"
  if (( code != 0 )); then
    tail -40 "$log_dir/$name.log"
    printf 'Failed: %s (exit %s). Logs: %s\n' "$name" "$code" "$log_dir" >&2
    exit "$code"
  fi
  printf 'Completed %s\n' "$name"
done
printf 'All 11 scripts completed. Logs: %s\n' "$log_dir"
