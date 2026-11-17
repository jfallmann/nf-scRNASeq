fixtures <- Sys.getenv("NFSCRNA_FIXTURES", unset = "/tmp/nfscrna-gex-fixtures")
dir.create(fixtures, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(Matrix)
    library(R.utils)
})

set.seed(42)
ngenes <- 208L
ncells <- 80L
marker.genes <- c("ACD", "ABCC5", "ANLN", "AHI1", "AGFG1", "CCND3", "ASF1B", "AURKB")
mt.genes <- c("MT-CO1", "MT-ND1", "MT-CYB", "MT-ATP6")
ribo.genes <- c("RPS6", "RPL13A", "RPS3", "RPL8")
bg.genes <- sprintf("G%03d", seq_len(ngenes - length(marker.genes) - length(mt.genes) - length(ribo.genes)))
genes <- c(marker.genes, mt.genes, ribo.genes, bg.genes)
stopifnot(length(genes) == ngenes, !anyDuplicated(genes))
de.genes <- sprintf("G%03d", 1:10)

gtf <- file.path(fixtures, "annotation.gtf")
gtf.lines <- sprintf(
    "chr1\tensembl\tgene\t%d\t%d\t.\t+\t.\tgene_id \"%s\"; gene_biotype \"%s\";",
    seq_len(ngenes) * 1000L, seq_len(ngenes) * 1000L + 500L, genes,
    ifelse(genes %in% mt.genes, "mitochondrial_gene", "protein_coding")
)
writeLines(gtf.lines, gtf)

libs <- data.frame(
    sample = c("S1_ctrl_1", "S2_ctrl_2", "S3_trt_1", "S4_trt_2"),
    condition = c("ctrl", "ctrl", "trt", "trt"),
    stringsAsFactors = FALSE
)

make_library <- function(i) {
    set.seed(100 + i)
    base <- matrix(rpois(ngenes * ncells, lambda = 3), nrow = ngenes, ncol = ncells,
        dimnames = list(genes, paste0("cell", seq_len(ncells), "-1")))
    base[mt.genes, ] <- rpois(length(mt.genes) * ncells, lambda = 1)
    if (libs$condition[i] == "trt") {
        base[de.genes, ] <- base[de.genes, , drop = FALSE] + rpois(length(de.genes) * ncells, lambda = 20)
    }
    Matrix(base, sparse = TRUE)
}

write_library <- function(i) {
    dir <- file.path(fixtures, libs$sample[i])
    unlink(dir, recursive = TRUE)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    m <- make_library(i)
    mtx <- file.path(dir, "matrix.mtx")
    writeMM(m, mtx)
    R.utils::gzip(mtx, destname = paste0(mtx, ".gz"), remove = TRUE, overwrite = TRUE)
    con <- gzfile(file.path(dir, "barcodes.tsv.gz"), "w")
    writeLines(colnames(m), con)
    close(con)
    con <- gzfile(file.path(dir, "features.tsv.gz"), "w")
    writeLines(paste(genes, genes, "Gene Expression", sep = "\t"), con)
    close(con)
    dir
}

lib.dirs <- vapply(seq_len(nrow(libs)), write_library, character(1))
stopifnot(all(file.exists(file.path(lib.dirs, "matrix.mtx.gz"))))
