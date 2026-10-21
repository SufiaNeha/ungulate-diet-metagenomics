#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:?Set TRIMMOMATIC_JAR}"

mkdir -p \
  results/trimmed \
  results/nonhost \
  results/logs

[[ -s config/adapters.fa ]] || {
  echo "Missing adapter FASTA" >&2
  exit 1
}

printf "SampleID\tNonhostPairs\n" \
  > results/logs/nonhost_pair_counts.tsv

while IFS=$'\t' read -r id r1 r2; do

  [[ "$id" == "SampleID" || -z "$id" || "$id" == \#* ]] &&
    continue

  [[ -s "$r1" && -s "$r2" ]] || {
    echo "Missing FASTQ: $id" >&2
    exit 1
  }

  java -jar "$TRIMMOMATIC_JAR" PE \
    -threads "$THREADS" \
    -phred33 \
    "$r1" "$r2" \
    "results/trimmed/${id}_R1_paired.fq.gz" \
    "results/trimmed/${id}_R1_unpaired.fq.gz" \
    "results/trimmed/${id}_R2_paired.fq.gz" \
    "results/trimmed/${id}_R2_unpaired.fq.gz" \
    ILLUMINACLIP:config/adapters.fa:2:30:10 \
    SLIDINGWINDOW:4:20 \
    MINLEN:90 \
    2> "results/logs/${id}.trimmomatic.log"

  # Keep a read pair only when BOTH mates are unmapped.
  bowtie2 \
    --very-sensitive \
    -x references/hosts \
    -p "$THREADS" \
    -1 "results/trimmed/${id}_R1_paired.fq.gz" \
    -2 "results/trimmed/${id}_R2_paired.fq.gz" \
    2> "results/logs/${id}.bowtie2.log" |
    samtools view \
      -@ "$THREADS" \
      -b -f 12 -F 0x900 - |
    samtools sort \
      -n -@ "$THREADS" \
      -o "results/nonhost/${id}.bam" -

  samtools fastq \
    -@ "$THREADS" \
    -n \
    -1 "results/nonhost/${id}_R1.fq.gz" \
    -2 "results/nonhost/${id}_R2.fq.gz" \
    -0 /dev/null \
    -s /dev/null \
    "results/nonhost/${id}.bam"

  n1=$(
    gzip -cd "results/nonhost/${id}_R1.fq.gz" |
      awk 'END {print NR/4}'
  )

  n2=$(
    gzip -cd "results/nonhost/${id}_R2.fq.gz" |
      awk 'END {print NR/4}'
  )

  [[ "$n1" == "$n2" ]] || {
    echo "Mate count mismatch: $id" >&2
    exit 1
  }

  printf "%s\t%s\n" "$id" "$n1" \
    >> results/logs/nonhost_pair_counts.tsv

done < config/samples.tsv
