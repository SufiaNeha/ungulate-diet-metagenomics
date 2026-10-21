#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
mkdir -p references results/logs

# Three FASTA paths, one per line:
# Aoudad: GCF_016772045.1
# Desert bighorn sheep: GCF_042477335.2
# Mule deer: GCA_020976825.1

mapfile -t genomes < <(
  grep -v '^#' config/host_genomes.txt |
  sed '/^[[:space:]]*$/d'
)

[[ ${#genomes[@]} -eq 3 ]] || {
  echo "Expected three host genomes" >&2
  exit 1
}

: > references/hosts_combined.fa

for i in "${!genomes[@]}"; do
  fasta="${genomes[$i]}"
  [[ -s "$fasta" ]] || exit 1

  awk -v prefix="host$((i+1))_" '
    /^>/ {
      sub(/^>/, ">" prefix)
      print
      next
    }
    {print}
  ' "$fasta" >> references/hosts_combined.fa
done

bowtie2-build \
  --threads "$THREADS" \
  references/hosts_combined.fa \
  references/hosts

# Build a plant-only Kraken2 database.
KRAKEN_DB="${KRAKEN_DB:-references/kraken_plant}"

if [[ "${BUILD_KRAKEN_DB:-0}" == "1" ]]; then
  kraken2-build --download-taxonomy --db "$KRAKEN_DB"
  kraken2-build --download-library plant --db "$KRAKEN_DB"
  kraken2-build --build \
    --threads "$THREADS" \
    --db "$KRAKEN_DB"

  bracken-build \
    -d "$KRAKEN_DB" \
    -t "$THREADS" \
    -k 35 \
    -l 90
fi
