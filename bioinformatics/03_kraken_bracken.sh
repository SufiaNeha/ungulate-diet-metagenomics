#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
KRAKEN_DB="${KRAKEN_DB:-references/kraken_plant}"

mkdir -p \
  results/kraken \
  results/bracken \
  results/logs

[[ -s "$KRAKEN_DB/hash.k2d" ]] || {
  echo "Kraken2 database missing" >&2
  exit 1
}

while IFS=$'\t' read -r id r1 r2; do

  [[ "$id" == "SampleID" || -z "$id" || "$id" == \#* ]] &&
    continue

  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired \
    --gzip-compressed \
    --report "results/kraken/${id}.report" \
    --output "results/kraken/${id}.kraken" \
    "results/nonhost/${id}_R1.fq.gz" \
    "results/nonhost/${id}_R2.fq.gz" \
    2> "results/logs/${id}.kraken.log"

  # Genus-level abundance
  bracken \
    -d "$KRAKEN_DB" \
    -i "results/kraken/${id}.report" \
    -o "results/bracken/${id}.G.tsv" \
    -r 90 \
    -l G \
    -t 0

  # Family-level abundance
  bracken \
    -d "$KRAKEN_DB" \
    -i "results/kraken/${id}.report" \
    -o "results/bracken/${id}.F.tsv" \
    -r 90 \
    -l F \
    -t 0

done < config/samples.tsv
