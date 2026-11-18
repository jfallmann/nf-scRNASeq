# nf-scRNASeq

Nextflow (DSL2) pipeline for single-cell RNA-seq libraries: mapping (CellRanger or STARsolo) or
reuse of precomputed CellRanger output, quality control, filtering, Seurat/SCE/AnnData object
creation, overview plots and replicate-level differential expression. It follows the CaTCHseq
pipeline up to gene-level analysis only: there is no CaTCH barcode counting/collapsing, so no
barcode parameters exist here.

## Requirements

- Nextflow (tested with 26.10.4) and a Java runtime.
- Docker (or Apptainer/Singularity through a Nextflow profile) to run the task containers.
- A mapper index and reference for the `STAR` route; `CellRanger` needs its own reference.

## Samplesheet

`--libraries` takes a CSV with these required columns:

| Column | Content |
| --- | --- |
| `SampleName` | sample name; repeat the row to spread a library over several lanes |
| `Condition` | condition of the sample (treatment, timepoint, genotype) |
| `Replicate` | positive integer, mandatory even for a single replicate |
| `Lane` | positive integer, orders the FASTQ files of a library |
| `LibraryType` | `GEX` (raw FASTQ) or `CellRanger` (precomputed CellRanger output directory given in `R1`) |
| `R1` | R1 FASTQ (`.fastq.gz`) or, for `CellRanger`, the CellRanger output directory |
| `R2` | R2 FASTQ (`.fastq.gz`); required for `PIPseq`, optional otherwise |
| `CellNumber` | expected cells: `NA`, a soft constraint (`7500`) or a hard constraint (`7500!`) |
| `Chemistry` | chemistry label, e.g. `10X`, `DropIn`, `SmartSeq`, `PIPseq` |

Raw GEX libraries:

```csv
SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry
S1,ctrl,1,1,GEX,/data/S1_L001_R1_001.fastq.gz,/data/S1_L001_R2_001.fastq.gz,7500,10X
S2,trt,1,1,GEX,/data/S2_L001_R1_001.fastq.gz,/data/S2_L001_R2_001.fastq.gz,7500,10X
```

Precomputed CellRanger output (the directory is copied, never modified in place):

```csv
SampleName,Condition,Replicate,Lane,LibraryType,R1,R2,CellNumber,Chemistry
S1,ctrl,1,1,CellRanger,/data/cellranger/S1,,NA,10X
S2,trt,1,1,CellRanger,/data/cellranger/S2,,NA,10X
```

Validate a samplesheet before running:

```bash
python3 docker/scripts/python/checkSampleSheet.py samplesheet.csv
```

## Chemistry presets

STARsolo options are defined once in `nextflow/conf/chemistry.config` under
`params.chemistries` and resolved by `nextflow/lib/Chemistry.groovy`. Each entry declares a
`mode`, a one-based `geometry` (`cbStart`, `cbLen`, `umiStart`, `umiLen`, plus `barcodeRead`
where the barcode sits in R2) and an ordered `starOptions` map. The geometry is authoritative:
the `--soloCBstart/--soloCBlen/--soloUMIstart/--soloUMIlen` options are taken from it.

Supported names in the `Chemistry` column: `10X` (generic, 16 bp CB + 12 bp UMI, the default for
any name containing `10X` that is not listed below), `10Xv2` (16 + 10), `10Xv3` (16 + 12),
`Droplet` (16 + 10), `ScaleBio` (29 + 8, barcode read 2), `Smart` (SmartSeq, no barcode geometry)
and `PIPseq` (24 + 12). `10x`, `droplet`, `scalebio`, `smart`, `smartseq` and `pipseq` are
accepted as aliases.

Override a preset in a `-c` config; assignments merge into the defaults, so only the changed
fields need to be given:

```groovy
params {
    chemistries {
        '10Xv2' {
            geometry {
                umiLen = 10
            }
        }
    }
}
```

The legacy overrides `star_10X`, `star_droplet`, `star_scalebio`, `star_smart` and `star_pipseq`
default to `null` and are applied when set: a string replaces the preset options; `star_pipseq`
additionally accepts a `{ cbLen, umiLen -> ... }` closure. `--pipseq_cb_length` and
`--pipseq_umi_length` change the converted PIPseq geometry and the UMI start follows the barcode
length; the raw preprocessing geometry (`--pipseq_offset`, `--pipseq_umi_offset`) is separate and
unchanged.

## Running

```bash
nextflow run nextflow/main.nf -with-docker \
    --libraries samplesheet.csv \
    --outputDir ./OUTPUT \
    --reportsDir ./REPORTS \
    --mapper CellRanger \
    --reference /data/ref.fa.gz \
    --annotation /data/genes.gtf.gz \
    --index /data/index \
    --baseline ctrl
```

`nextflow run nextflow/main.nf --help` prints the authoritative parameter list with the live
defaults. The analysis defaults are:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--output_seurat` | `true` | publish Seurat objects |
| `--output_sce` | `false` | build and publish `SingleCellExperiment` objects |
| `--output_anndata` | `false` | publish AnnData `.h5ad` objects |
| `--run_sct` | `true` | run SCTransform and its reductions/clustering |
| `--run_pathways` | `false` | decoupleR/progeny pathway inference (needs network access, not exercised by the baseline tests) |
| `--de_rds` | `null` | run only the DE step on precomputed Seurat `.rds`/`.rds.gz` files (comma separated); no mapping |

Disabling an export flag only stops publication: the filtered Seurat object is still built
internally and still feeds the overview plots and the DE step.

## PIPseq

PIPseq libraries need `--mapper STAR`; the preprocessing runs in
`ghcr.io/obenauflab/pipspeak:latest` (`--pipseq_container`). The bundled configuration
`nextflow/conf/pipseq/pipspeak.yaml` (`--pipseq_config`) is the default and resolves its barcode
lists relative to the YAML file; those lists are staged automatically, cached deeply and the input
files are left untouched. For VBCF `24CB12UMI` barcodes use the defaults
`--pipseq_cb_length 24 --pipseq_umi_length 12 --pipseq_offset 8 --pipseq_umi_offset 10`.

## Output

```
OUTPUT/
├── RAWMAPPING/
│   ├── CellRanger/            # per-library CellRanger output
│   └── STAR/                  # per-library STARsolo output
├── QC/
│   ├── FASTQC/                # --withQC
│   └── MULTI/                 # --withQC
├── PIPseq/                    # PIPseq logs, QC tables and converted <lib>_R1/R2.fastq.gz
├── SCE/
│   ├── filtered/              # <lib>_nf-scRNASeq_filtered_seurat_sce.rds.gz
│   │                          # <lib>_nf-scRNASeq_filtered_sce.rds.gz
│   └── raw/                   # <lib>_nf-scRNASeq_unfiltered_seurat_sce.rds.gz
│                              # <lib>_nf-scRNASeq_unfiltered_sce.rds.gz
├── AnnData/
│   ├── filtered/              # <lib>_nf-scRNASeq_filtered.h5ad
│   └── raw/                   # <lib>_nf-scRNASeq_unfiltered.h5ad
├── Plots/
│   └── Overview/<lib>/        # <lib>_nf-scRNASeq_overview.pdf
│                              # <lib>_nf-scRNASeq_overview_QC_summary.tsv
│                              # <lib>_nf-scRNASeq_overview_cell_metadata.tsv
│                              # <lib>_nf-scRNASeq_overview_QC_*.pdf
└── DE/
    ├── GENES/                 # allSamples_nf-scRNASeq_DE_summary.rds.gz
    │                          # allSamples_nf-scRNASeq_DE_skipped.tsv
    │                          # DE_DESEQ2_<contrast>_GENES_all.tsv.gz
    │                          # DE_DESEQ2_<contrast>_GENES_significant.tsv.gz
    └── Plots/                 # DE_DESEQ2_<contrast>_GENES_Volcano_*.pdf
```

`<lib>` is `SampleName_Condition_Replicate` joined by underscores, with underscores inside each
field replaced by dashes. With `--filter true` the object names carry the `nf-scRNASeq.prefiltered`
infix instead of `nf-scRNASeq`.

## Differential expression

DE is computed on replicate-level pseudobulk counts (DESeq2), one contrast per non-baseline
condition against `--baseline`. A contrast is skipped, not fabricated, when the baseline condition
is absent, when only one condition is present, or when a condition has fewer than two replicates;
each skip is recorded with its reason in `DE/GENES/allSamples_nf-scRNASeq_DE_skipped.tsv`. No
cluster-level DE is produced.

## Containers

All analysis tasks run in `ghcr.io/obenauflab/catchseq:latest`, overridable with `--r_container`.
The image supplies the shared R/Bioconductor dependencies while the pipeline's own R scripts are
staged into each task from `--scriptDirR`; the shared helpers (`singlecell_utils.R`,
`object_exports.R`) are adapted from the committed CaTCHseq versions. Mapping tasks use
`nfcore/cellranger:7.1.0` and `ghcr.io/obenauflab/catchseq:latest` for STAR; the container mapping
lives in `nextflow/conf/docker.config`. Building the image locally is optional — the published
image is enough to run the pipeline.

## Tests

```bash
python3 tests/test_workflow.py          # end-to-end workflow on synthetic fixtures
Rscript tests/test_gex_processing.R     # R script suite (needs the catchseq R environment)
python3 tests/test_pipseq.py            # PIPseq helper and samplesheet tests
```

The fixtures are generated deterministically by `tests/make_fixtures.R`; set `NFSCRNA_FIXTURES` to
reuse an existing fixture directory instead of the default `/tmp/nfscrna-gex-fixtures`.
