# nf-neoag

A [Nextflow](https://www.nextflow.io/) (DSL2) pipeline for neoantigen prediction,
covering:

- neoantigens from SNVs** (`neo-snv`)
- neoantigens from Indels** (`neo-indel`)
- *(optional)* merge of the two consensus outputs (`neo-merge`)

The analysis is a self-contained, dependency-free **Python 3 package** (`neoag`,
in [`src/neoag`](src/neoag)) that reimplements the original method. The three
`neo-*` commands can also be run standalone, without Nextflow.

> **Attribution.** The neoantigen-prediction **method** is the work of
> **Giuseppe Rospo**, a former group member of Genomics and Targeted Therapies Lab led by Prof. Alberto Bardelli  — all scientific credit goes to him.
> This repository is an independent Python 3 reimplementation of that method,
> wrapped as a Nextflow pipeline for portability and reproducibility.
> See [`CITATIONS.md`](CITATIONS.md).

## Authors & contributions

- **Giuseppe Rospo** — original neoantigen-prediction method and reference
  implementation (the scientific core). Please cite Germano G *et al.*,
  *Nature* (2017); see [`CITATIONS.md`](CITATIONS.md).
- **Vittorio Battaglieri** — Python 3 reimplementation (`neoag` package),
  Nextflow port, containerization, multi-genome configuration (hg38/hg19/mm10),
  CI, and packaging for reproducible, shareable use.


## What it does

For each sample, starting from **already-called variants** (not from FASTQ/BAM):

```
SNV table  (.snp.str-bias.Idea.OK)  ──▶  NEO_SNV   ──▶  <sample>.neoantigens.results.consensus
Indel table (.cosmic.Indel.altX)    ──▶  NEO_INDEL ──▶  <sample>.neoantigens.results.consensus
                                               │
                                (--merge)      └──▶  <sample>.neoantigens.consensus.all[.expanded]
```

Peptides (8–11-mers) are generated around each variant, scored with **netMHC-4.0**,
and collapsed to a per-variant consensus.

## Requirements

- Nextflow ≥ 22.10
- Docker or Singularity
- **netMHC-4.0** — bring your own (license-restricted, not bundled)
- **Reference files** per genome — bring your own (large, not bundled). See
  [`docs/references.md`](docs/references.md).

Supported genomes: **hg38, hg19, mm10**.
(hg19 currently has no `.yar` reference, so for hg19 only the **SNV** step is available.)

Reference files already present on this machine can be wired up with the
provided (gitignored) local config:

```bash
nextflow run . -profile docker -c conf/refs.local.config \
  --input samplesheet.csv --netmhc /opt/netMHC-4.0
```

## Quick smoke test (no netMHC / no real refs)

Wiring test using the two bundled example inputs and `-stub-run`:

```bash
nextflow run . -profile test,docker -stub-run
```

## Real run

1. Build/pull the container (build context is the repo root, so the `neoag`
   package can be installed into the image):
   ```bash
   docker build -t nf-neoag:0.1.0 -f containers/Dockerfile .
   ```
   (or set `--container` / `process.container` to a published image)

2. Prepare a samplesheet (see [`assets/samplesheet_example.csv`](assets/samplesheet_example.csv)):

   ```csv
   sample,genome,snv_table,indel_table,hla
   patientA,hg38,/data/A.snp.str-bias.Idea.OK,/data/A.cosmic.Indel.alt5,"HLA-A0201,HLA-B0702"
   ```

   - Leave `snv_table` empty to run only the indel step (and vice-versa).
   - `hla` is a netMHC allele list; **quote it** if it contains commas.

3. Run:

   ```bash
   nextflow run . -profile docker \
     --input samplesheet.csv \
     --outdir results \
     --netmhc /opt/netMHC-4.0 \
     --genomes.hg38.cdna_tab /data/refs/hg38/refFlat_mRNA.hg38.noalt.parsed.tab \
     --genomes.hg38.yar      /data/refs/hg38/refFlat_mRNA.hg38.yar \
     --merge
   ```

   (Or set the reference paths once in `conf/genomes.config`.)

## Parameters

| param             | default            | description                                             |
|-------------------|--------------------|---------------------------------------------------------|
| `--input`         | —                  | samplesheet CSV (required)                              |
| `--outdir`        | `results`          | output directory                                        |
| `--netmhc`        | —                  | netMHC-4.0 install **directory** (required)             |
| `--netmhc_bin`    | `netMHC`           | executable name inside `--netmhc`                       |
| `--merge`         | `false`            | also merge SNV+Indel consensus per sample               |
| `--genomes.<b>.cdna_tab` | (in config) | `refFlat_mRNA.<b>.noalt.parsed.tab` (SNV step)          |
| `--genomes.<b>.yar`      | (in config) | `refFlat_mRNA.<b>.yar` (Indel step)                     |

## Running the analysis without Nextflow

The `neoag` package installs three standalone commands. Point `NETMHC` (or
`--netmhc`) at your netMHC-4.0 executable/directory:

```bash
pip install .

export NETMHC=/opt/netMHC-4.0/netMHC
neo-snv   sample.snp.str-bias.Idea.OK  refFlat_mRNA.hg38.noalt.parsed.tab  sample  "HLA-A0201,HLA-B0702"
neo-indel sample.cosmic.Indel.alt5     refFlat_mRNA.hg38.yar               sample  "HLA-A0201,HLA-B0702"
neo-merge sample.neoantigens.results.consensus  sample.neoantigens.results.consensus  sample
```

Each writes `sample.neoantigens.results.consensus`; `neo-merge` writes
`sample.neoantigens.consensus.all[.expanded]`.


## Layout

```
main.nf                    workflow
modules/local/*.nf         NEO_SNV, NEO_INDEL, MERGE_SNV_FS
src/neoag/                 Python 3 analysis package (neo-snv/neo-indel/neo-merge)
pyproject.toml             package metadata + console entrypoints
conf/                      base / genomes / test configs
containers/Dockerfile      Python 3 + tcsh runtime (installs the neoag package)
assets/                    samplesheet schema + example
test/                      example inputs + stub fixtures
docs/references.md         how to obtain the reference files
```
