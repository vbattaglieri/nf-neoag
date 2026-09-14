# Reference files

Two reference files are required **per genome build**. They are large and are
**not** shipped with this pipeline. You supply them via `conf/genomes.config`
(or `--genomes.<build>.<key>` on the command line, or a custom `-c` config).

| key        | file                                      | used by            | source (original SOP paths)                                  |
|------------|-------------------------------------------|--------------------|--------------------------------------------------------------|
| `cdna_tab` | `refFlat_mRNA.<build>.noalt.parsed.tab`   | step 10 (neo/SNV)  | `/mnt/Unito-Analysis/big/Giuseppe/IDEA_package/refFlat_mRNA.<build>.noalt.parsed.tab` |
| `yar`      | `refFlat_mRNA.<build>.yar`                | step 11 (neo/Indel)| `/scratch/reference/UCSC_<build>/ANNOYAR/refFlat_mRNA.<build>.yar` |

Supported builds: **hg38, hg19, mm10**.

Availability of the files currently at hand:

| build | `cdna_tab` (SNV) | `yar` (Indel) |
|-------|:----------------:|:-------------:|
| hg38  | yes              | yes           |
| hg19  | yes              | **missing**   |
| mm10  | yes              | yes           |

hg19 therefore supports the SNV step only until an hg19 `.yar` is produced/obtained.

## How to obtain them

The scripts that originally **built** these two files are not available, so they
cannot be regenerated from scratch here. Recommended path:

1. **Copy the existing files** off the lab cluster (paths above) to a stable
   location you control.
2. Point the pipeline at them, e.g.:

   ```bash
   nextflow run . -profile docker \
     --input samplesheet.csv \
     --netmhc /opt/netMHC-4.0 \
     --genomes.hg38.cdna_tab /data/refs/hg38/refFlat_mRNA.hg38.noalt.parsed.tab \
     --genomes.hg38.yar      /data/refs/hg38/refFlat_mRNA.hg38.yar
   ```

3. **(Recommended for publication)** Once you have confirmed the files, deposit
   them on **Zenodo** (gives a citable DOI). Nextflow can stage them directly
   from the resulting URLs — just set the URL as the value in `genomes.config`.
   The repository then stays small and anyone can reproduce the run.

> These files must never be committed to git (see `.gitignore`).

## netMHC-4.0

netMHC is license-restricted and is therefore **not** included in the container.
Install it yourself, then pass its directory with `--netmhc /path/to/netMHC-4.0`.
The executable name inside that directory is assumed to be `netMHC`
(override with `--netmhc_bin`).
