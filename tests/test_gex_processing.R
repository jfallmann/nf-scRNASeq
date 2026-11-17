root <- Sys.getenv("NFSCRNA_ROOT", unset = getwd())
script.dir <- file.path(root, "docker", "scripts", "R")
marker.file <- file.path(root, "docker", "data", "R", "stagemarkers_xue2020.rds")
stopifnot(dir.exists(script.dir), file.exists(marker.file))

fixtures <- Sys.getenv("NFSCRNA_FIXTURES", unset = "/tmp/nfscrna-gex-fixtures")
outdir <- "/tmp/nfscrna-gex-out"
logdir <- "/tmp/nfscrna-gex-logs"
dir.create(fixtures, showWarnings = FALSE, recursive = TRUE)
unlink(outdir, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(logdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
    library(SingleCellExperiment)
})
options(Seurat.object.assay.version = "v5")

source(file.path(root, "tests", "make_fixtures.R"))

run_cli <- function(script, args, tag) {
    log <- file.path(logdir, paste0(tag, ".log"))
    wd <- file.path(outdir, tag)
    dir.create(wd, showWarnings = FALSE, recursive = TRUE)
    old <- getwd()
    setwd(wd)
    on.exit(setwd(old), add = TRUE)
    status <- system2("Rscript",
        c("--vanilla", file.path(script.dir, script), args),
        stdout = log, stderr = paste0(log, ".err")
    )
    if (!identical(as.integer(status), 0L)) {
        cat(readLines(log, warn = FALSE), sep = "\n")
        cat(readLines(paste0(log, ".err"), warn = FALSE), sep = "\n")
        stop(paste("CLI failed:", script, tag))
    }
    invisible(status)
}

preprocess_args <- function(i, out, output_seurat, output_sce, output_anndata, run_sct) {
    c(
        "--sample", libs$sample[i],
        "--data10X", lib.dirs[i],
        "--annotation", gtf,
        "--baseCond", "ctrl",
        "--max_mt", "100",
        "--min_features", "10",
        "--hvg_cutoff", "0.1",
        "--out", out,
        "--marker", marker.file,
        "--libpath", paste0(script.dir, "/"),
        "--output_seurat", output_seurat,
        "--output_sce", output_sce,
        "--output_anndata", output_anndata,
        "--run_sct", run_sct
    )
}

out1 <- file.path(outdir, "lib1_all")
run_cli("preprocessData.R", preprocess_args(1, out1, "true", "true", "true", "true"), "preprocess_lib1")

out2 <- file.path(outdir, "lib2_minimal")
run_cli("preprocessData.R", preprocess_args(2, out2, "false", "false", "false", "false"), "preprocess_lib2")

out3 <- file.path(outdir, "lib3_default")
run_cli("preprocessData.R", preprocess_args(3, out3, "true", "false", "false", "true"), "preprocess_lib3")

out4 <- file.path(outdir, "lib4_default")
run_cli("preprocessData.R", preprocess_args(4, out4, "true", "false", "false", "true"), "preprocess_lib4")

cat("preprocessData.R CLI runs OK\n")

pdf_pages <- function(path) {
    raw <- readBin(path, "raw", file.size(path))
    hits <- grepRaw("/Type /Page", raw, all = TRUE, fixed = TRUE)
    sum(raw[hits + 11L] != charToRaw("s"))
}

expected1 <- paste0(out1, c(
    "_unfiltered_seurat_sce.rds.gz", "_unfiltered_sce.rds.gz", "_unfiltered.h5ad",
    "_filtered_seurat_sce.rds.gz", "_filtered_sce.rds.gz"
))
stopifnot(all(file.exists(expected1)))

expected2 <- paste0(out2, c(
    "_unfiltered_seurat_sce.rds.gz", "_unfiltered_sce.rds.gz", "_unfiltered.h5ad",
    "_filtered_sce.rds.gz"
))
stopifnot(!any(file.exists(expected2)))
stopifnot(file.exists(paste0(out2, "_filtered_seurat_sce.rds.gz")))
cat("output flag contract OK\n")

seu <- readRDS(paste0(out1, "_filtered_seurat_sce.rds.gz"))
stopifnot(nrow(seu) == ngenes, ncol(seu) == ncells)
stopifnot(setequal(rownames(seu), genes))
meta.cols <- c(
    "Sample", "Condition", "Replicate", "percent.mt", "percent.ribo", "Phase", "CellStage",
    "pca_cluster", "sct_cluster", "is.damaged", "is.low_yield", "Category"
)
stopifnot(all(meta.cols %in% colnames(seu@meta.data)))
stopifnot(identical(as.character(unique(seu$Sample)), "S1"))
stopifnot(identical(as.character(unique(seu$Condition)), "ctrl"))
stopifnot(identical(as.character(unique(seu$Replicate)), "ctrl_1"))
stopifnot(all(c("RNA", "SCT") %in% names(seu@assays)))
stopifnot(all(c("pca", "pca_sct", "umap_pca", "umap_pca_sct") %in% names(seu@reductions)))
rna.layers <- Layers(seu[["RNA"]])
stopifnot(any(grepl("^counts", rna.layers)), any(grepl("^data", rna.layers)))
counts.layer <- LayerData(seu[["RNA"]], layer = grep("^counts", rna.layers, value = TRUE)[1])
stopifnot(is(counts.layer, "sparseMatrix"), identical(dim(counts.layer), c(ngenes, ncells)))
stopifnot(all(colnames(counts.layer) == colnames(seu)))
stopifnot(all(seu$percent.mt >= 0), all(seu$percent.ribo >= 0))
cat("filtered Seurat readback OK\n")

sce <- readRDS(paste0(out1, "_filtered_sce.rds.gz"))
stopifnot(is(sce, "SingleCellExperiment"))
stopifnot(all(c("counts", "logcounts") %in% assayNames(sce)))
stopifnot(identical(dim(sce), c(ngenes, ncells)))
stopifnot(all(c("GeneID", "Symbol", "Biotype", "Locus", "is.mitochondrial", "is.ribosomal") %in% colnames(rowData(sce))))
stopifnot(all(c("Sample", "Condition", "Replicate", "Cluster", "CellStage", "Category") %in% colnames(colData(sce))))
stopifnot(is(assay(sce, "counts"), "sparseMatrix"))
cat("filtered SCE readback OK\n")

ad <- anndataR::read_h5ad(paste0(out1, "_unfiltered.h5ad"))
stopifnot(identical(as.integer(dim(ad)), c(ncells, ngenes)))
seu.unf <- readRDS(paste0(out1, "_unfiltered_seurat_sce.rds.gz"))
ad.data <- LayerData(JoinLayers(seu.unf[["RNA"]]), layer = "data")
stopifnot(identical(rownames(ad$X), colnames(ad.data)), identical(colnames(ad$X), rownames(ad.data)))
stopifnot(isTRUE(all.equal(as.matrix(t(as.matrix(ad$X))), as.matrix(ad.data))))
stopifnot(identical(rownames(ad$obs), colnames(seu.unf)), identical(rownames(ad$var), rownames(seu.unf)))
cat("unfiltered AnnData readback OK\n")

seu2 <- readRDS(paste0(out2, "_filtered_seurat_sce.rds.gz"))
stopifnot(nrow(seu2) == ngenes, ncol(seu2) == ncells)
stopifnot(all(c("pca", "umap_pca") %in% names(seu2@reductions)))
stopifnot(!("pca_sct" %in% names(seu2@reductions)), !("SCT" %in% names(seu2@assays)))
cat("minimal flag run readback OK\n")

overview.out <- file.path(outdir, "lib1_overview")
run_cli("create_overview_plots.R", c(
    "--sce", paste0(out1, "_filtered_seurat_sce.rds.gz"),
    "--baseCond", "ctrl",
    "--out", overview.out,
    "--format", "pdf",
    "--width", "25",
    "--height", "10",
    "--libpath", paste0(script.dir, "/")
), "overview_lib1")

overview.pdf <- paste0(overview.out, ".pdf")
stopifnot(file.exists(overview.pdf), file.size(overview.pdf) > 0)
stopifnot(pdf_pages(overview.pdf) >= 2)
qc.summary <- read.delim(paste0(overview.out, "_QC_summary.tsv"), check.names = FALSE)
stopifnot(nrow(qc.summary) == 1L, all(c("Sample", "Condition", "Replicate", "n_cells") %in% colnames(qc.summary)))
cell.meta <- read.delim(paste0(overview.out, "_cell_metadata.tsv"), check.names = FALSE)
stopifnot(nrow(cell.meta) == ncells)
cat("overview plots OK (", pdf_pages(overview.pdf), "pages )\n")

sce.paths <- function(outs) {
    paste(vapply(outs, function(o) paste0(o, "_filtered_seurat_sce.rds.gz"), character(1)), collapse = ",")
}

de.out <- file.path(outdir, "de_main")
run_cli("identify_de_genes.R", c(
    "--sce", sce.paths(c(out1, out2, out3, out4)),
    "--baseCond", "ctrl",
    "--run_pathways", "false",
    "--pcut", "0.1",
    "--fcut", "1",
    "--out", de.out,
    "--libpath", paste0(script.dir, "/")
), "de_main")

de.dir <- file.path(outdir, "de_main")
all.file <- file.path(de.dir, "DE_DESEQ2_trt-vs-ctrl_GENES_all.tsv.gz")
sig.file <- file.path(de.dir, "DE_DESEQ2_trt-vs-ctrl_GENES_significant.tsv.gz")
volcano.file <- file.path(de.dir, "DE_DESEQ2_trt-vs-ctrl_GENES_Volcano_p0.1_lfc1.pdf")
stopifnot(file.exists(all.file), file.exists(sig.file), file.exists(volcano.file))
stopifnot(pdf_pages(volcano.file) >= 1)

deres <- read.delim(gzfile(all.file), check.names = FALSE)
stopifnot(nrow(deres) == ngenes, all(c("Gene", "log2FoldChange", "log2FoldChange_MLE", "p.adj") %in% colnames(deres)))
stopifnot(all(deres$log2FoldChange_MLE[deres$Gene %in% de.genes] > 0))
stopifnot(mean(deres$log2FoldChange_MLE[deres$Gene %in% de.genes]) > mean(deres$log2FoldChange_MLE[!deres$Gene %in% de.genes]))

de.summary <- readRDS(paste0(de.out, "_DE_summary.rds.gz"))
stopifnot(nrow(de.summary$pseudobulk_samples) == 4L)
stopifnot(!("cluster" %in% colnames(de.summary$pseudobulk_samples)))
stopifnot(setequal(de.summary$pseudobulk_samples$Replicate, c("ctrl_1", "ctrl_2", "trt_1", "trt_2")))
stopifnot(setequal(de.summary$pseudobulk_samples$Condition, c("ctrl", "trt")))
stopifnot(identical(de.summary$contrasts, "trt-vs-ctrl"))

de.skipped <- read.delim(paste0(de.out, "_DE_skipped.tsv"), check.names = FALSE)
stopifnot(nrow(de.skipped) == 0L)
pathway.skipped <- read.delim(paste0(de.out, "_DE_pathways_skipped.tsv"), check.names = FALSE)
stopifnot(nrow(pathway.skipped) == 1L, grepl("disabled", pathway.skipped$reason))
comparisons <- readRDS(paste0(de.out, "_DE_pseudobulk_list.rds.gz"))
stopifnot(identical(names(comparisons), "trt-vs-ctrl"))
cat("DE replicate pseudobulk OK\n")

skip.single <- file.path(outdir, "de_single_condition")
run_cli("identify_de_genes.R", c(
    "--sce", sce.paths(c(out1, out2)),
    "--baseCond", "ctrl",
    "--run_pathways", "false",
    "--out", skip.single,
    "--libpath", paste0(script.dir, "/")
), "de_single_condition")
sk.single <- read.delim(paste0(skip.single, "_DE_skipped.tsv"), check.names = FALSE)
stopifnot(nrow(sk.single) == 1L, grepl("single_condition", sk.single$contrast))
stopifnot(!file.exists(file.path(outdir, "de_single_condition", "DE_DESEQ2_trt-vs-ctrl_GENES_all.tsv.gz")))

skip.reps <- file.path(outdir, "de_insufficient_reps")
run_cli("identify_de_genes.R", c(
    "--sce", sce.paths(c(out1, out3)),
    "--baseCond", "ctrl",
    "--run_pathways", "false",
    "--out", skip.reps,
    "--libpath", paste0(script.dir, "/")
), "de_insufficient_reps")
sk.reps <- read.delim(paste0(skip.reps, "_DE_skipped.tsv"), check.names = FALSE)
stopifnot(nrow(sk.reps) == 1L, identical(sk.reps$contrast, "trt-vs-ctrl"), grepl("Fewer than 2 replicates", sk.reps$reason))
stopifnot(!file.exists(file.path(outdir, "de_insufficient_reps", "DE_DESEQ2_trt-vs-ctrl_GENES_all.tsv.gz")))

skip.base <- file.path(outdir, "de_missing_baseline")
run_cli("identify_de_genes.R", c(
    "--sce", sce.paths(c(out1, out2, out3, out4)),
    "--baseCond", "nope",
    "--run_pathways", "false",
    "--out", skip.base,
    "--libpath", paste0(script.dir, "/")
), "de_missing_baseline")
sk.base <- read.delim(paste0(skip.base, "_DE_skipped.tsv"), check.names = FALSE)
stopifnot(nrow(sk.base) == 1L, grepl("missing_baseline", sk.base$contrast))
stopifnot(!file.exists(file.path(outdir, "de_missing_baseline", "DE_DESEQ2_trt-vs-ctrl_GENES_all.tsv.gz")))
cat("DE skip cases OK\n")

cat("ALL TESTS PASSED\n")
