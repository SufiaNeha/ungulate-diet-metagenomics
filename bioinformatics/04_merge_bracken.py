#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

outdir = Path("results/matrices")
outdir.mkdir(parents=True, exist_ok=True)

with open("config/samples.tsv", newline="") as f:
    samples = [
        row["SampleID"]
        for row in csv.DictReader(f, delimiter="\t")
        if row["SampleID"]
        and not row["SampleID"].startswith("#")
    ]

if not samples or len(samples) != len(set(samples)):
    raise ValueError("Sample IDs must be present and unique")

for rank, label in [("G", "genus"), ("F", "family")]:

    counts = defaultdict(dict)
    names = {}

    for sample in samples:
        filename = Path(
            f"results/bracken/{sample}.{rank}.tsv"
        )

        with filename.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                taxid = row["taxonomy_id"]
                abundance = float(row["new_est_reads"])

                counts[taxid][sample] = (
                    counts[taxid].get(sample, 0.0)
                    + abundance
                )

                names[taxid] = row["name"]

    output = outdir / f"{label}_counts_by_taxid.csv"

    with output.open("w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow(
            ["taxid", "Taxon"] + samples
        )

        for taxid in sorted(counts, key=int):
            writer.writerow(
                [taxid, names[taxid]]
                + [
                    counts[taxid].get(sample, 0)
                    for sample in samples
                ]
            )

    print(
        f"{label}: {len(counts)} taxa, "
        f"{len(samples)} samples"
    )
