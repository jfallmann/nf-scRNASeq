.libPaths("")

library(optparse)

option_list <- list(
    make_option(c("--sample"),
        type = "character", default = NULL,
        help = "sample name"
    ),
    make_option(c("--baseCond"),
        type = "character", default = NULL,
        help = "Base condition for DE Analysis"
    ),
    make_option(c("--data10X"),
        type = "character", default = NULL,
        help = "path to the 10X data"
    ),
    make_option(c("--annotation"),
        type = "character", default = NULL,
        help = "path to the matching annotation file in GTF format"
    ),
    make_option(c("--max_mt"),
        type = "numeric", default = 10,
        help = "maximum percent of mitochondrial reads in a valid cell"
    ),
    make_option(c("--min_features"),
        type = "numeric", default = 500,
        help = "minimum number of detected features in a valid cell"
    ),
    make_option(c("--hvg_cutoff"),
        type = "numeric", default = 0.1,
        help = "cutoff value to call percentage of high variable genes (must be between 0 and 1)"
    ),
    make_option(c("--output_seurat"), type = "character", default = "true",
        help = "save unfiltered Seurat output (true/false); filtered Seurat is retained internally"
    ),
    make_option(c("--output_sce"), type = "character", default = "false",
        help = "process and save SingleCellExperiment output (true/false)"
    ),
    make_option(c("--output_anndata"), type = "character", default = "false",
        help = "save unfiltered AnnData output (true/false)"
    ),
    make_option(c("--run_sct"), type = "character", default = "true",
        help = "run SCTransform and SCT reductions/clustering (true/false)"
    ),
    make_option(c("--out"),
        type = "character", default = NULL,
        help = "path to the output file"
    ),
    make_option(c("--marker"),
        type = "character", default = "/tools/data/R/stagemarkers_xue2020.rds",
        help = "path to cellcycle marker file in RDS format"
    ),
    make_option(c("--libpath"),
        type = "character", default = "/tools/scripts/R/",
        help = "path to R libs, trailing slash is required"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(paste0("CLI: ", opt))

if (is.null(opt$sample) || is.null(opt$out) || is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("Not enough parameters")
}

#### Source Functions ####
source(paste0(opt$libpath, "singlecell_utils.R"))
source(paste0(opt$libpath, "object_exports.R"))
for (name in c("output_seurat", "output_sce", "output_anndata", "run_sct")) {
    opt[[name]] <- parse_object_flag(opt[[name]], name)
}

############################

library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(Seurat)
# options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
if (opt$run_sct) library(sctransform)
library(R.filesets)

### Load all experiments as layer and converting to Seuratv5 Object
sl <- create_SCEs(opt$sample, opt$data10X, opt$annotation)

sce <- if (opt$output_sce) sl$sce else NULL
seurat_sce <- sl$seurat_sce

print(paste("Seurat object contains", ncol(seurat_sce), "cells and", nrow(seurat_sce), "genes."))

#### Remove cells without any counts from both objects ####
nonzero <- colnames(sl$sce)[colSums(counts(sl$sce)) > 0]
rm(sl) # clean up

if (length(nonzero) == 0) {
    stop("No cell with a non-zero count remains after loading the 10X data")
}
if (!is.null(sce)) {
    sce <- sce[, nonzero]
}
seurat_sce <- subset(seurat_sce, cells = nonzero)
print(paste("Keeping", length(nonzero), "cells with non-zero counts."))

#### Annotate Samples, Conditions and Replicates ####
if (opt$output_sce) {
    sce$Condition <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
    sce$Replicate <- as.factor(unlist(lapply(sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
    sce$Sample <- as.factor(unlist(lapply(sce$Sample, function(x) unlist(str_split(x, "_"))[1])))
}

seurat_sce$Sample.orig <- as.factor(seurat_sce$Sample)
seurat_sce$Condition <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[2])))
seurat_sce$Replicate <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) paste(unlist(str_split(x, "_"))[2], unlist(str_split(x, "_"))[3], sep = "_"))))
seurat_sce$Sample <- as.factor(unlist(lapply(seurat_sce$Sample, function(x) unlist(str_split(x, "_"))[1])))

#### calculate percentage of MT reads SEURAT ####
seurat_sce <- PercentageFeatureSet(seurat_sce, pattern = "^MT-|^mt-", col.name = "percent.mt")
seurat_sce <- PercentageFeatureSet(seurat_sce, pattern = "^RP[SL]|^rp[sl]|^Rp[sl]", col.name = "percent.ribo")
pdf(file = paste0(opt$out, "_QC_Violin_MT_content.pdf"), width = 30, height = 10)
VlnPlot(seurat_sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_MT_content.pdf"), width = 30, height = 10)
FeatureScatter(seurat_sce, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
dev.off()
#
pdf(file = paste0(opt$out, "_QC_Scatter_Feature_content.pdf"), width = 30, height = 10)
FeatureScatter(seurat_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
dev.off()

#### normalize counts and calculate percentage of MT reads SCE ####
if (opt$output_sce) {
    print("Normalize sce ...")

    is.mitochondrial <- grepl(
        pattern = "^MT-|^mt-", # Needed e.g. for mouse 10x data with cellranger prebuilt index
        x = rownames(sce),
        ignore.case = FALSE,
        perl = TRUE
    )
    is.ribosomal <- grepl(
        pattern = "^RP[SL]|^rp[sl]|^Rp[sl]",
        x = rownames(sce),
        ignore.case = FALSE,
        perl = TRUE
    )

    sce <- sce %>%
        scater::logNormCounts() %>%
        scater::addPerCellQC(subsets = list(MT = is.mitochondrial, RB = is.ribosomal)) %>%
        scater::addPerFeatureQC()
}

### Normalize and scale counts Seurat ####
print("Normalize Seurat ...")

seurat_sce <- normalize_Seurat(seurat_sce)

#### annotate low yield and damaged cells ####
if (opt$output_sce) {
    print("Annotate sce ...")
    rowData(sce)["is.mitochondrial"] <- is.mitochondrial
    rowData(sce)["is.ribosomal"] <- is.ribosomal
    colData(sce)["is.damaged"] <- colData(sce)[, "subsets_MT_percent"] > opt$max_mt
    colData(sce)["is.low_yield"] <- colData(sce)[, "detected"] < opt$min_features
}

print("Annotate seurat ...")
seurat_sce[["is.damaged"]] <- seurat_sce@meta.data %>%
    pull(percent.mt) %>%
    {
        case_when(. > opt$max_mt ~ TRUE, .default = FALSE)
    }
seurat_sce[["is.low_yield"]] <- seurat_sce@meta.data %>%
    pull(nFeature_RNA) %>%
    {
        case_when(. < opt$min_features ~ TRUE, .default = FALSE)
    }

#### Categorize the cells ####
cell.categories <- c("Good", "Damaged", "Few features", "Damaged AND few features")

if (opt$output_sce) {
    print("Categorize sce ...")
    colData(sce)["Category"] <- colData(sce) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
        dplyr::select("tmp") %>%
        purrr::map(.f = ~ cell.categories[.x + 1]) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
        dplyr::select(Class)
}

print("Categorize seurat ...")

seurat_sce@meta.data$Category <- seurat_sce@meta.data %>%
    dplyr::mutate(tmp = as.integer(is.damaged) * 1 + as.integer(is.low_yield) * 2) %>%
    dplyr::select("tmp") %>%
    purrr::map(.f = ~ cell.categories[.x + 1]) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Class = factor(tmp, levels = cell.categories)) %>%
    dplyr::pull(Class)

### Assign cell stage and categories
print(paste0("Loading cell stage markers from ", opt$marker, sep = ""))
markerfile <- loadRDS(opt$marker)

#### Joining seurat layers for downstream annotation ####
print("Joining seurat layers ...")
seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "data", new = "data")
seurat_sce <- join_Seurat(seurat_sce, assay = "RNA", layers = "counts", new = "counts")

#### Attempting to assign cell phase via Seurat ####
print("Phase annotation seurat ...")

s.genes <- markerfile$S
g2m.genes <- markerfile$G2M

phase.assigned <- TRUE
tryCatch(
    withCallingHandlers(
        {
            seurat_sce <- CellCycleScoring(seurat_sce, assay = "RNA", slot = "data", s.features = str_to_upper(s.genes), g2m.features = str_to_upper(g2m.genes), set.ident = FALSE)
        },
        warning = function(w) {
            if (grepl("Could not find enough features in the object", w$message)) {
                print("WARNING COUGHT:  Could not find enough features in the object. Will try to match case")
                seurat_sce <- CellCycleScoring(seurat_sce, assay = "RNA", slot = "data", s.features = str_to_title(s.genes), g2m.features = str_to_title(g2m.genes), set.ident = FALSE)
            } else {
                message(w$message)
            }
        }
    ),
    error = function(e) {
        print(paste("ERROR COUGHT:  ", e, " WILL SET SEURAT CELL PHASE TO UNASSIGNED"))
        phase.assigned <<- FALSE
    }
)

if (!phase.assigned) {
    seurat_sce@meta.data$Phase <- NA_character_
    seurat_sce@meta.data$S.Score <- NA_real_
    seurat_sce@meta.data$G2M.Score <- NA_real_
}

#### Attempting to assign cell stage to SCEs ####
if (opt$output_sce) {
    print("Cellstage sce ...")
    sce <- sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage")
}
print("Cellstage seurat ...")
seurat_sce <- seurat_sce %>% assignCategoryByMarker(markers = markerfile, col.name = "CellStage")

### Split SCE again
seurat_sce <- split_Seurat(seurat_sce, by = seurat_sce$Sample.orig)

#### Save unfiltered sce and seurat_sce objects ####
print("Save unfiltered ...")
if (opt$output_sce) {
    saveRDS(sce, file = paste0(opt$out, "_unfiltered_sce.rds.gz"), compress = "gzip")
}
if (opt$output_seurat) {
    saveRDS(seurat_sce, file = paste0(opt$out, "_unfiltered_seurat_sce.rds.gz"), compress = "gzip")
}
if (opt$output_anndata) {
    write_anndata(seurat_sce, paste0(opt$out, "_unfiltered.h5ad"))
}

### Filter for MT content and min reads
print("Filter seurat ...")
# seurat_sce <- subset(seurat_sce, subset = nFeature_RNA > opt$min_features & percent.mt < opt$max_mt)
seurat_sce <- subset(seurat_sce, subset = is.low_yield == FALSE & is.damaged == FALSE)
if (opt$output_sce) {
    print("Filter sce ...")
    sce <- sce[, sce$is.low_yield == FALSE & sce$is.damaged == FALSE]
}
if (ncol(seurat_sce) == 0) {
    stop("No cell passed the filtering for mitochondrial content and minimum number of features")
}
print(paste0("Keeping ", ncol(seurat_sce), " Cells from Seurat object."))


#### Run PCA and UMAP ####
if (opt$output_sce) {
    print("Identify the top variable genes...")
    gene.var <- modelGeneVar(sce)
    hvg <- getTopHVGs(stats = gene.var, prop = opt$hvg_cutoff)

    print("Run the PCA ...")
    sce <- runPCA(sce, subset_row = hvg)
    set.seed(42)

    print("Run tSNE and UMAP analyses ...")
    sce <- runTSNE(sce, dimred = "PCA")
    sce <- runUMAP(sce, dimred = "PCA")

    print("Clustering ...")
    g <- buildSNNGraph(sce, use.dimred = "PCA")
    cluster <- igraph::cluster_walktrap(g)$membership
    colData(sce)["Cluster"] <- factor(cluster)
}


### SCtransform counts
if (opt$run_sct) {
    print("SCtransform Seurat ...")
    seurat_sce <- sctransform_Seurat(seurat_sce)
}
### Run initial dim reduction for norm
print("PCA Seurat ...")
seurat_sce <- reduceDims_Seurat(seurat_sce)
### Run initial dim reduction for SCT
if (opt$run_sct) {
    seurat_sce <- reduceDims_Seurat(seurat_sce, assay = "SCT", reduction.name = "pca_sct")
}
### Cluster
print("Clustering Seurat ...")
seurat_sce <- cluster_Seurat(seurat_sce, assay = "RNA", reduction = "pca", cluster.name = "pca_cluster")
if (opt$run_sct) {
    seurat_sce <- cluster_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", cluster.name = "sct_cluster")
}
## Run UMAPs
print("UMAP Seurat ...")
seurat_sce <- umap_Seurat(seurat_sce, assay = "RNA", reduction = "pca", dims = 1:30, reduction.name = "umap_pca", n.neighbors = 30L, min.dist = 0.1, spread = 5)
if (opt$run_sct) {
    seurat_sce <- umap_Seurat(seurat_sce, assay = "SCT", reduction = "pca_sct", dims = 1:30, reduction.name = "umap_pca_sct", n.neighbors = 30L, min.dist = 0.1, spread = 5)
}

#### Save final objects ####
print("Final Save ...")

if (opt$output_sce) {
    saveRDS(sce, file = paste0(opt$out, "_filtered_sce.rds.gz"), compress = "gzip")
}
saveRDS(seurat_sce, file = paste0(opt$out, "_filtered_seurat_sce.rds.gz"), compress = "gzip")
if (opt$output_anndata) {
    write_anndata(seurat_sce, paste0(opt$out, "_filtered.h5ad"))
}
