#!/usr/bin/env nextflow

/************************************************************************
===========================
nf-scRNASeq
===========================

Authors
- Joerg Fallmann (joerg.fallmann@imp.ac.at)
************************************************************************/

nextflow.enable.dsl=2

//normalize potentially Boolean/null string params to avoid concatenation issues
def normalizeStr(value){
    return (value instanceof Boolean || value == null) ? '' : value.toString()
}

def canonicalLibraryType(value){
    def text = value?.toString()?.trim()
    if (text?.toLowerCase() == 'gex') return 'GEX'
    if (text?.toLowerCase() == 'cellranger') return 'CellRanger'
    if (text?.toLowerCase() == 'catchseq') error "LibraryType 'CaTCHseq' is not supported by nf-scRNASeq, use GEX"
    error "LibraryType '${text}' is not supported by nf-scRNASeq, use GEX or CellRanger"
}

def normalizeChemistry(value){
    def text = value?.toString()?.trim()
    def lower = text?.toLowerCase()
    if (lower == 'pipseq') return 'PIPseq'
    if (lower == '10x') return '10X'
    if (lower == 'droplet') return 'Droplet'
    if (lower == 'scalebio') return 'ScaleBio'
    if (lower == 'smart' || lower == 'smartseq') return 'Smart'
    return text
}

def sq(value){
    def text = value.toString()
    if (text.contains("'")) error "path '${text}' must not contain a single quote"
    return "'" + text + "'"
}

def validate_pipseq_geometry(){
    def cbLen = params.pipseqCbLength.toString().toInteger()
    def umiLen = params.pipseqUmiLength.toString().toInteger()
    def offset = params.pipseqOffset.toString().toInteger()
    def umiOffset = params.pipseqUmiOffset.toString().toInteger()
    def clip3p = (params.pipseqClip3p != null) ? params.pipseqClip3p.toString().toInteger() : 0
    if (cbLen <= 0) error "--pipseq_cb_length must be positive, got ${cbLen}"
    if (umiLen < 1 || umiLen > 16) error "--pipseq_umi_length must be between 1 and 16, got ${umiLen}"
    if (offset < 0) error "--pipseq_offset must not be negative, got ${offset}"
    if (umiOffset < 0) error "--pipseq_umi_offset must not be negative, got ${umiOffset}"
    if (clip3p < 0) error "--pipseq_clip3p must not be negative, got ${clip3p}"
}

def resolvePipspeakConfig(value){
    def configFile = file(value, checkIfExists: true)
    if (!configFile.isFile()) error "PIPseq config '${configFile}' is not a regular file"
    def options = new org.yaml.snakeyaml.LoaderOptions()
    options.setAllowDuplicateKeys(false)
    def raw
    try {
        raw = new org.yaml.snakeyaml.Yaml(new org.yaml.snakeyaml.constructor.SafeConstructor(options)).load(configFile.text)
    } catch (Throwable failure) {
        error "PIPseq config '${configFile}' is not valid YAML: ${failure.message}"
    }
    if (!(raw instanceof Map)) error "PIPseq config '${configFile}' does not contain a YAML mapping"
    def barcodes = raw['barcodes']
    if (!(barcodes instanceof Map) || barcodes.isEmpty()) error "PIPseq config '${configFile}' has no 'barcodes' mapping"
    def names = []
    def lists = []
    barcodes.each { name, reference ->
        if (!(name instanceof CharSequence) || name.toString().trim().isEmpty()) {
            error "PIPseq config '${configFile}' has a barcode set without a name"
        }
        if (!(reference instanceof CharSequence) || reference.toString().trim().isEmpty()) {
            error "PIPseq config '${configFile}' barcode list '${name}' has no path"
        }
        def relative = reference.toString()
        if (relative.startsWith('/')) {
            error "PIPseq config '${configFile}' barcode list '${name}' path '${relative}' must be relative to the config file"
        }
        def listFile = configFile.parent.resolve(relative).normalize()
        if (!listFile.isFile()) {
            error "PIPseq config '${configFile}' barcode list '${name}' path '${relative}' does not resolve to a file (${listFile})"
        }
        names.add(name.toString())
        lists.add(listFile)
    }
    return [configFile, names, lists]
}

def objectFlag(value, name) {
    def text = value?.toString()?.toLowerCase()
    if (!(text in ['true', 'false'])) {
        error "--${name} must be true or false"
    }
    return text == 'true'
}

def flagOn(value, fallback) {
    def text = value?.toString()?.toLowerCase()
    return text == null ? fallback : text == 'true'
}

def mappingPublishPath(filename, subtree){
    def parts = filename.tokenize('/')
    def head = parts[0]
    def suffixes = ['_filtered_feature_bc_matrix', '_raw_feature_bc_matrix', '_filtered_barcodes.tsv.gz', '_raw_barcodes.tsv.gz', '_mapped.sam.gz', '_mapped.bam.bai', '_mapped.bam', '_mapping.log', '_SJ.out.tab', '_R1_unmapped.fastq.gz', '_R2_unmapped.fastq.gz']
    def lib = head
    def match = suffixes.find { head.endsWith(it) }
    if (match) {
        lib = head.substring(0, head.length() - match.length())
    }
    def rest = (parts.size() > 1) ? parts[1..-1].join('/') : head
    return "${subtree}/${lib}/${rest}"
}

//Params from CL
//most params are used directly as params.<name> (see conf/default-params.config for defaults).
//the following are aliases for CLI flags whose names differ from their internal default names.
params.absDir = workflow.launchDir
params.qcparams = normalizeStr(params.fastqc_params ?: params.qcParams)
params.mapperbin = params.mapper ?: params.mapperBin
params.runqc = params.withQC ?: params.runQC
params.crparams = normalizeStr(params.cellranger_params ?: params.crParams)
params.starparams = normalizeStr(params.star_params ?: params.starMappingParams)
params.idxparams = normalizeStr(params.idx_params ?: params.idxParams)
params.mapindex = params.index ?: params.mapIndex
params.mapref = params.reference ?: params.mapRef
params.mapanno = params.annotation ?: params.mapAnno
params.markerfile = params.markers ?: params.markerFile
params.maxmtpercent = params.max_mt_percent ?: params.maxMtPercent
params.mindetectedfeatures = params.min_detected_features ?: params.minDetectedFeatures
params.hvgcutoff = params.hvg_cutoff ?: params.hvgCutoff
params.pvalcutoff = params.pval_cutoff ?: params.pvalCutoff
params.lfccutoff = params.lfc_cutoff ?: params.lfcCutoff
params.refName = params.baseline ?: 'Day0'
params.rContainer = params.r_container ?: 'ghcr.io/obenauflab/catchseq:latest'
params.output_seurat = objectFlag(params.containsKey('output_seurat') && params.output_seurat != null ? params.output_seurat : true, 'output_seurat')
params.output_sce = objectFlag(params.containsKey('output_sce') && params.output_sce != null ? params.output_sce : false, 'output_sce')
params.output_anndata = objectFlag(params.containsKey('output_anndata') && params.output_anndata != null ? params.output_anndata : false, 'output_anndata')
params.run_sct = objectFlag(params.containsKey('run_sct') && params.run_sct != null ? params.run_sct : true, 'run_sct')
params.run_pathways = objectFlag(params.containsKey('run_pathways') && params.run_pathways != null ? params.run_pathways : false, 'run_pathways')


def helpMessage() {
    log.info"""
    ======================================================================
      nf-scRNASeq

      Singlecell RNA-Seq pipeline


      Version: ${workflow.manifest.version}
      Contact: Joerg Fallmann (joerg.fallmann@imp.ac.at)
    ======================================================================

      Usage:
      nextflow run nf-scRNASeq/nextflow/main.nf --libraries <SAMPLESHEET>

      Mandatory arguments:
        --libraries             CSV file with the following columns:
                                    SampleName      name of the sample (can appear in multiple lines, in case
                                                    the library was sequenced in several runs)
                                    Condition       condition for this sample (e.g. timepoint, KO, treatment)
                                    Replicate       replicate (even if a single replicate is present, this column cannot be missing or be empty)
                                    Lane            sequencing lane, used to order the FASTQ files of a library
                                    LibraryType     either GEX (raw FASTQ) or CellRanger (precomputed CellRanger output directory in R1)
                                    R1              path to the R1 read, or to a precomputed CellRanger output directory
                                    R2              path to the R2 read (if available)
                                    CellNumber      number of expected cells (NA if not available, number for soft constraint, number! for hard constraint)
                                    Chemistry       chemistry used (e.g. 10X, DropIn, SmartSeq, PIPseq), this may not set adequate parameters for mappers automatically, make sure you check them accordingly

            Optional arguments:
                --outputDir             specifies the output directory (default: ${params.outputDir})
                --reportsDir            specifies the reports directory.(default: ${params.reportsDir})
                --scriptDirR            specifies the path to the R scripts directory, it is staged into every task (default: ${params.scriptDirR})
                --scriptDirPy           specifies the path to the Python scripts directory, it is staged into every task (default: ${params.scriptDirPy})
                --binDir                specifies the path to the binary directory, do not change if running with docker, otherwise set to path on nf-scRNASeq git repo.(default: /usr/local/bin/ which is valid for docker instance; set to '' for instances not running docker)
                --mapper                Which mapper to run (default: CellRanger, optional: STAR)
                --index                 Path to mapper index directory (default: ${params.mapindex}, NEEDS TO BE SET ALSO TO CREATE NEW INDEX, new index will be stored at given path)
                --reference             Path to reference fasta.gz
                --annotation            Path to annotation gtf.gz
                --whitelist             Path to barcode whitelist
                --withQC                Boolean, run FastQC and MultiQC (default: ${params.runqc})
                --fastqc_params         Optional parameters for FASTQC
                --cellranger_params     Optional parameters for CellRanger (default: ${params.crparams})
                --star_params           Optional parameters for STAR mapping
                --idx_params            Optional parameters for STAR index generation
                --filter                Postprocess filtered counts (default: ${params.filter})
                --minReads              Minimum reads to keep cell (default: ${params.minReads})
                --max_mt_percent        Maximum mitochondrial read percentage (default: ${params.maxmtpercent})
                --min_detected_features Minimum detected features cutoff (default: ${params.mindetectedfeatures})
                --hvg_cutoff            Highly variable genes cutoff (default: ${params.hvgcutoff})
                --output_seurat         Publish Seurat objects (true/false, default: ${params.output_seurat})
                --output_sce            Process/publish SingleCellExperiment objects (true/false, default: ${params.output_sce})
                --output_anndata        Publish AnnData .h5ad objects (true/false, default: ${params.output_anndata})
                --run_sct               Run SCTransform and SCT reductions/clustering (true/false, default: ${params.run_sct})
                --pval_cutoff           P-value cutoff for DE analysis (default: ${params.pvalcutoff})
                --lfc_cutoff            Log-fold-change cutoff for DE analysis (default: ${params.lfccutoff})
                --run_pathways          Run decoupleR/progeny pathway inference, needs network access (true/false, default: ${params.run_pathways})
                --de_rds                Run only the DE analysis on these precomputed Seurat .rds/.rds.gz files (comma separated), no mapping is run (default: ${params.de_rds})
                --organism              Identifier for organism (choice: ["Human", "Mouse"], default: ${params.organism})
                --baseline              Name of reference day/condition (default: ${params.refName})
                --marker                RDS file of cellcycle markers (default: ${params.markerfile}, NamedList with gene names for each stage [G1S, S, G2M, M, MG1, G0], S and G2M are needed)
                --pipseq_config         PIPseq barcode configuration YAML, its barcode lists are resolved relative to the YAML file (default: ${params.pipseqConfig})
                --pipseq_cb_length      Length of the concatenated PIPseq barcode (default: ${params.pipseqCbLength})
                --pipseq_umi_length     Length of the PIPseq UMI (default: ${params.pipseqUmiLength})
                --pipseq_offset         Barcode offset from the start of R1 (default: ${params.pipseqOffset})
                --pipseq_umi_offset     UMI offset from the end of the barcode (default: ${params.pipseqUmiOffset})
                --pipseq_exact          Use exact barcode matching in pipspeak (true/false, default: ${params.pipseqExact})
                --pipseq_clip3p         STAR --clip3pNbases for PIPseq libraries (default: ${params.pipseqClip3p})
                --pipseq_container      Container image for the PIPseq preprocessing (default: ${params.pipseqContainer})
                --r_container           Container image for the R analysis processes (default: ${params.rContainer})
                --stopOnWarnings        Stop pipeline when warnings occur (default: ${params.stopOnWarnings})
                --help                  print this help message

    """.stripIndent()
}


// ---------------------------------------------------------------------
// Pipeline Channels and Processes
// ---------------------------------------------------------------------

// For more information about syntax, please refer to the nextflow documentation at https://www.nextflow.io/docs/latest/index.html

/************************************************************************
                SANITY CHECK SAMPLESHEET
************************************************************************/

process check_samplesheet {
    tag "$samplesheet"

    input:
    path samplesheet
    path checker

    output:
    path 'Valid_*.csv', emit: csv
    path 'sanity_check', emit: check

    script: 
    """
    python3 ${checker} \
        ${samplesheet} \
    && cp -L ${samplesheet} Valid_${samplesheet}
    """
}


/************************************************************************
                CONCATENATE LANES
************************************************************************/

process concat_lanes {

    tag { sampleName instanceof List ? sampleName.join(',') : sampleName }

    input:
    tuple val(sampleName), path(r1, stageAs: 'inputs/r1??/*'), path(r2, stageAs: 'inputs/r2??/*'), val(cells_expected), val(chemistry)

    output:
    tuple val(sampleName), path("${sampleName}_R1.fastq.gz",includeInputs:false), path("${sampleName}_R2.fastq.gz",includeInputs:false), val(cells_expected), val(chemistry)

    script:
    """
    cat \$(ls inputs/r1??/* | sort) > "${sampleName}_R1.fastq.gz"
    cat \$(ls inputs/r2??/* | sort) > "${sampleName}_R2.fastq.gz"
    """
}


/************************************************************************
                STEP 0.0: Prepare PIPseq libraries
************************************************************************/

process preprocessPipseq {

    cache 'deep'
    tag "${meta.id}"

    publishDir "${params.outputDir}/PIPseq", mode: 'copy'

    input:
    tuple val(meta), path(r1, stageAs: 'inputs/r1??/*'), path(r2, stageAs: 'inputs/r2??/*'), path(configYaml, stageAs: 'config/pipspeak.yaml'), val(barcodeNames), path(barcodeFiles, stageAs: 'barcodes/list??/*'), path(helper)

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: converted
    path("pipspeak/*_log.yaml"), emit: logs, optional: true
    path("pipspeak/*.tsv"), emit: qc, optional: true
    path("pipspeak/*_whitelist.txt"), emit: whitelist, optional: true
    path("pipspeak/*_cutadapt.log"), emit: cutadapt_log, optional: true

    script:
    def exactFlag = objectFlag(params.pipseqExact, 'pipseq_exact') ? '--exact' : ''
    def r1Args = (r1 instanceof List ? r1 : [r1]).collect { sq(it) }.join(' ')
    def r2Args = (r2 instanceof List ? r2 : [r2]).collect { sq(it) }.join(' ')
    def stagedLists = (barcodeFiles instanceof List ? barcodeFiles : [barcodeFiles])
    if (stagedLists.size() != barcodeNames.size()) {
        error "PIPseq config lists ${stagedLists.size()} barcode files for ${barcodeNames.size()} barcode sets"
    }
    def barcodeArgs = [barcodeNames, stagedLists].transpose().collect { name, staged -> "--barcode-file ${sq(name)} ${sq(staged)}" }.join(' ')
    """
    python3 ${sq(helper)} \
        --r1 ${r1Args} \
        --r2 ${r2Args} \
        --prefix ${sq(meta.id)} \
        --config ${sq(configYaml)} \
        ${barcodeArgs} \
        --cb-len ${params.pipseqCbLength} \
        --umi-len ${params.pipseqUmiLength} \
        --offset ${params.pipseqOffset} \
        --umi-offset ${params.pipseqUmiOffset} \
        --threads ${task.cpus} \
        --library-type ${sq(meta.librarytype)} \
        ${exactFlag}
    """
}


/************************************************************************
                STEP 0: Run read QC
************************************************************************/

process qc_raw{
   
    cache 'lenient'
    tag "${sampleName}"

    publishDir "${params.outputDir}/QC/FASTQC", mode: 'copy'
    
    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    path("*.zip"), emit: zip
    path("*.html"), emit: html

    script:
    if (params.binDir){
        fqc = params.binDir+"fastqc"
    } else{
        fqc = "fastqc"
    }
    """
    ${fqc} --quiet -t ${task.cpus} $params.qcparams --noextract -f fastq $read1 $read2 &&
    for fqc in *_fastqc.{zip,html}
    do
        mv "\$fqc" "${sampleName}_\$fqc"
    done
    """
}

process mqc{
    cache 'lenient'

    publishDir "${params.outputDir}/QC/MULTI", mode: 'copy'

    input:
    path(fastqcs)

    output:
    path "*.zip", includeInputs:false, emit: mqc
    path "*.html", includeInputs:false, emit: html

    script:
    if (params.binDir != ''){
        mqc = params.binDir+"multiqc"
    } else{
        fqc = "multiqc"
    }
    """
    touch $fastqcs; export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; ${mqc} -f -k json -z -o \${PWD} .
    """
}

/************************************************************************
                STEP 1: Run CellRanger count
************************************************************************/

process Cellranger_idx{

    cache 'lenient'
    label 'big_mem'

    publishDir "${params.outputDir}/RAWMAPPING/CellRanger/LOGS", mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)       "${file(filename).getName()}"
        else if (filename.indexOf(".idx") > 0)     "${params.mapindex}.idx"
        else                                       "${params.mapindex}"
    }

    input:
    path genome
    path anno

    output:
    path "${IDX}", emit: idx
    path "${IDX}*", emit: idx_extra
    path "star.idx", emit: idxlink
    path "*.out", emit: idxlog

    script:
    gen =  file(genome).getName()
    an  = file(anno).getName()
    IDX = file(params.mapindex).getName()
    taskmem = task.memory.toGiga()

    """
    zcat ${gen} > tmp.fa \
        && zcat ${an} > tmp_anno \
        && cellranger mkref ${params.idxparams}\
        --genome=${IDX} \
        --nthreads ${task.cpus} \
        --memgb ${taskmem} \
        --fasta tmp.fa \
        --genes tmp_anno \
        && rm -f tmp.fa tmp_anno \
        && ln -s ${IDX}/star star.idx
    """
}


process runCellrangerCount{

    cache 'lenient'
    tag "${sampleName}"

    publishDir "${params.outputDir}", mode: 'copy',
    saveAs: {filename -> mappingPublishPath(filename, 'RAWMAPPING/CellRanger')}

    input:
        tuple val(sampleName), path('inputs/'), path('inputs/'), val(cells_expected), val(chemistry), path(index)

    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_raw
        tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
        tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    
    script:
    // Check chemistry specific settings
    chemistry = chemistry.unique()[0]
    cells_expected = cells_expected.unique()[0]
    if (!chemistry.contains("10X")){
        log.error("Running CellRanger on chemistry different than 10X is not supported, please check your settings and sample sheet.")

    }
    crparams = params.crparams
    if (cells_expected != "NA"){
        if (cells_expected.contains("!")){
            crparams = crparams + ' --force-cells ' + cells_expected.replaceAll('!', '')
        }else{
            crparams = crparams + ' --expect-cells ' + cells_expected 
        }
    }
    taskmem = task.memory.toGiga()

    """
    # Find all reads, sort them by name to ensure that the paired files are on the consecutive lines,
    # and then create symlinks with proper names (SampleName_S1_R1_xxx.fastq.gz)
    IDX=1
    BKP=\${IFS}
    IFS=\$'\\n'
    for LINE in \$(find inputs/ -regextype posix-extended -regex ".*R[12][\\._].*fastq.gz" -exec readlink -f {} \\; | sort | paste - -);
    do
        SUFFIX=\$(printf "%03d" \${IDX})

        mkdir -p tomap
        R1=\$(echo \${LINE} | cut -f1|sed 's|*/||g')
        NEW_NAME=${sampleName}_S1_R1_\${SUFFIX}.fastq.gz
        ln -sf \${R1} tomap/\${NEW_NAME}

        R2=\$(echo \${LINE} | cut -f2|sed 's|*/||g')
        NEW_NAME=${sampleName}_S1_R2_\${SUFFIX}.fastq.gz
        ln -sf \${R2} tomap/\${NEW_NAME}

        IDX=\$((IDX + 1))
    done

    IFS=\${BKP}

    cellranger count \
        --disable-ui \
        ${crparams} \
        --jobmode local \
        --localcores ${task.cpus} \
        --localmem ${taskmem} \
        --transcriptome ${index} \
        --id ${sampleName} \
        --fastqs tomap

    mv ${sampleName} rundir
    mv rundir/outs ${sampleName}
    ln -fs ${sampleName}/filtered_feature_bc_matrix ${sampleName}_filtered_feature_bc_matrix
    ln -fs ${sampleName}/raw_feature_bc_matrix ${sampleName}_raw_feature_bc_matrix
    """
}


process useCellrangerData{

    cache 'lenient'
    tag "${sampleName}"

    // No publishDir here: the source data the user points R1 at is, by definition,
    // already where they want it. Republishing risked writing back onto that same
    // path (via publishDir's default overwrite) and destroying the original.

    input:
        tuple val(sampleName), path("cr_data")

    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_from_precomputed_raw
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix"), emit: cell_data_from_precomputed_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix"), emit: cell_data_from_precomputed_raw

    script:
        """
        # Materialize real, independent copies of only what we need. Never write
        # into "cr_data" itself: it may be (or contain) a symlink straight back to
        # the user's original, externally-owned CellRanger output directory.
        mkdir -p ${sampleName}
        cp -rL cr_data/filtered_feature_bc_matrix ${sampleName}/filtered_feature_bc_matrix
        cp -rL cr_data/raw_feature_bc_matrix ${sampleName}/raw_feature_bc_matrix
        ln -fs ${sampleName}/filtered_feature_bc_matrix ${sampleName}_filtered_feature_bc_matrix
        ln -fs ${sampleName}/raw_feature_bc_matrix ${sampleName}_raw_feature_bc_matrix
        zcat ${sampleName}_raw_feature_bc_matrix/barcodes.tsv.gz > ${sampleName}_raw_feature_bc_matrix/barcodes.tsv
        """
}

/************************************************************************
                STEP 1 (Optional): Run STARsolo count
************************************************************************/

process star_idx{

    cache 'lenient'
    label 'big_mem'

    publishDir "${params.outputDir}/RAWMAPPING/STAR/LOGS", mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)       "${file(filename).getName()}"
        else if (filename.indexOf(".idx") > 0)     "${params.mapindex}.idx"
        else                                       "${params.mapindex}"
    }

    input:
    path genome
    path anno

    output:
    path "${IDX}", emit: idx
    path "${IDX}.idx", emit: idxlink
    path "*.out", emit: idxlog

    script:
    gen =  genome.getName()
    an  = anno.getName()
    IDX = file(gen).getSimpleName()+'_idx'
    taskmem      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    zcat ${gen} > tmp.fa && zcat ${an} > tmp_anno && mkdir -p ${IDX} && STAR ${params.idxparams} --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir ${IDX} --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno ${taskmem} && mv -f ${IDX}/*.out ${IDX}.Log.out && ln -s ${IDX} ${IDX}.idx
    """
}

process create_dummy_whitelist{
    cache 'lenient'
    label 'big_mem'

    output:
    path("dummy.txt"), emit: dummy

    script:
    """
    touch dummy.txt
    """
}

process star_mapping{
    cache 'lenient'
    label 'big_mem'
    tag "${sampleName}"

    publishDir "${params.outputDir}", mode: 'copy',
    saveAs: {filename -> mappingPublishPath(filename, 'RAWMAPPING/STAR')}

    input:
    tuple val(sampleName), path(r1), path(r2), val(cells_expected), val(chemistry), path(idx), path(whitelist)
    
    output:
    path "${sampleName}", emit: name
    tuple val(sampleName), path("${sampleName}"), emit: out
    tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
    tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    tuple val(sampleName), path("${sampleName}_filtered_barcodes.tsv.gz"), emit: cell_ids_filtered
    tuple val(sampleName), path("${sampleName}_raw_barcodes.tsv.gz"), emit: cell_ids_raw
    tuple val(sampleName), path("*_mapped.sam.gz"), emit: sam
    tuple val(sampleName), path("*.bam"), emit: bam
    tuple val(sampleName), path("*.bai"), emit: bai
    tuple val(sampleName), path("*.log"), emit: logs
    tuple val(sampleName), path("*Log*.out"), emit: xtralogs
    tuple val(sampleName), path("*.tab"), emit: sjtab
    tuple val(sampleName), path("*_unmapped.fastq.gz", includeInputs:false), emit: unmapped

    script:
    idxdir = idx.toRealPath()
    extraparams = ''
    starparams = params.starparams
    
    chemistry = chemistry.unique()[0]
    cells_expected = cells_expected.unique()[0]
    
    // Check whitelist
    if( (whitelist.size() > 0 ) && (chemistry != 'ScaleBio') && (chemistry != 'PIPseq')){
        starparams = starparams + " --soloCBwhitelist ${whitelist}"
    }else{
        starparams = starparams + " --soloCBwhitelist None"
    }
    // Check chemistry specific settings
    if (chemistry == "Droplet"){
        log.info("Running StarSolo on unspecified Droplet chemistry, please ensure your STARsolo parameters fit the protocol, please check and adapt default settings in mappers.config file in the conf directory of the nextflow subdirectory of this pipeline. There is no automatic sanity check!")
    } else if (chemistry == 'ScaleBio'){
        starparams = params.starMappingParamsScale
    } else if (chemistry == "Smart"){
        log.info("Running StarSolo on chemistry different than Droplet based, please ensure your STARsolo parameters fit the protocol, please check and adapt default settings in mappers.config file in the conf directory of the nextflow subdirectory of this pipeline. There is no automatic sanity check!")
    }
    extraparams = Chemistry.starArgs(Chemistry.resolve(chemistry, params))
    // Check expected cell count
    if (cells_expected != "NA"){
        cellsexpCR = cells_expected.replaceAll('!', '') + ' 0.99 10'
        cellsexpED = cells_expected.replaceAll('!', '') + ' 0.99 10 45000 90000 500 0.01 20000 0.01 10000'
        starparams = starparams.replaceAll('CellRanger2.2', 'CellRanger2.2 '+cellsexpCR).replaceAll('EmptyDrops_CR', 'EmptyDrops_CR '+cellsexpED)
        extraparams = extraparams.replaceAll('CellRanger2.2', 'CellRanger2.2 '+cellsexpCR).replaceAll('EmptyDrops_CR', 'EmptyDrops_CR '+cellsexpED)
    } else if (chemistry == 'PIPseq'){
        extraparams = extraparams.replaceAll('EmptyDrops_CR', 'EmptyDrops_CR 5000 0.99 10 45000 90000 500 0.01 20000 0.001 10000')
    }
    if (chemistry == 'PIPseq'){
        starparams = starparams.replaceAll(/--soloMultiMappers\s+\S+/, '')
    }
    
    // Build params
    starparams = starparams + ' ' + extraparams

    if (chemistry == 'PIPseq' && params.pipseqClip3p != null && params.pipseqClip3p.toString().toInteger() > 0){
        starparams = starparams + " --clip3pNbases ${params.pipseqClip3p}"
    }

    // Check read order
    if ( starparams.contains('--soloBarcodeMate' )){
        if ( starparams.contains('--soloBarcodeMate 1' )){
            read1 = r2
            read2 = r1
            starparams = starparams.replaceAll('--soloBarcodeMate 1', '' )
        }else if ( starparams.contains('--soloBarcodeMate 2' )){
            read1 = r1
            read2 = r2
            starparams = starparams.replaceAll('--soloBarcodeMate 2', '' )
        }else{
            error('specified --soloBarcodeMate with unknown read')
        }
    }else{
        read1 = r2
        read2 = r1
    }

    //Get specific Feature counts, we use the first feature in --soloFeatures as default
    if ( chemistry == 'PIPseq' ){
        sfeature = 'GeneFull'
    }else if ( starparams.contains('--soloFeatures') ){
        sfeature = starparams.split('--soloFeatures')[1].split(' ')[1]
    }else{
        sfeature = 'Gene'
    }

    of = sampleName+'.Aligned.sortedByCoord.out.bam'
    gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")
    gb = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.bam")

    """
    STAR ${starparams} --runThreadN ${task.cpus} --genomeDir ${idxdir} --readFilesCommand zcat --readFilesIn <(cat ${read1}) <(cat ${read2}) --outFileNamePrefix ${sampleName}. --outReadsUnmapped Fastx &&samtools view -h ${of} | gzip > ${gf} && touch ${sampleName}.Unmapped.out.mate1 ${sampleName}.Unmapped.out.mate2 && cat ${sampleName}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R1_unmapped.fastq.gz && cat ${sampleName}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R2_unmapped.fastq.gz && mv *Log.out ${sampleName}_mapping.log

    mv ${of} ${gb}
    samtools index ${gb}
    mv ${sampleName}.Solo.out ${sampleName}
    gzip -f ${sampleName}/*/filtered/*
    gzip -f ${sampleName}/*/raw/*
    ln -s ${sampleName}/${sfeature}/filtered ${sampleName}_filtered_feature_bc_matrix
    ln -s ${sampleName}/${sfeature}/raw ${sampleName}_raw_feature_bc_matrix
    ln -s ${sampleName}/${sfeature}/filtered/barcodes.tsv.gz ${sampleName}_filtered_barcodes.tsv.gz
    ln -s ${sampleName}/${sfeature}/raw/barcodes.tsv.gz ${sampleName}_raw_barcodes.tsv.gz
    """

}


/************************************************************************
                    STEP 2: Generate SingleCellExperiment/Seurat objects
************************************************************************/

process preprocessSingleCellData{

    cache 'lenient'
    tag "Preprocess_SingleCellData_${sampleName}"

    publishDir "${params.outputDir}", mode: 'copy',
    saveAs: {filename ->
        def name = file(filename).getName()
        def lib = name.replaceFirst(/_(nf-scRNASeq)(\.prefiltered)?_.*$/, '')
        if (name.endsWith('_filtered_seurat_sce.rds.gz')) {
            flagOn(params.output_seurat, true) ? "SCE/filtered/${name}" : null
        } else if (name.endsWith('_unfiltered_seurat_sce.rds.gz')) {
            flagOn(params.output_seurat, true) ? "SCE/raw/${name}" : null
        } else if (name.endsWith('_filtered_sce.rds.gz')) {
            "SCE/filtered/${name}"
        } else if (name.endsWith('_unfiltered_sce.rds.gz')) {
            "SCE/raw/${name}"
        } else if (name.endsWith('_filtered.h5ad')) {
            "AnnData/filtered/${name}"
        } else if (name.endsWith('_unfiltered.h5ad')) {
            "AnnData/raw/${name}"
        } else if (name.endsWith('.pdf') || name.endsWith('.tsv') || name.endsWith('.tsv.gz')) {
            "Plots/Overview/${lib}/${name}"
        } else {
            null
        }
    }

    input:
        tuple val(sampleName), path(featureMatrix), path(gtf), path(marker), path(rscripts)

    output:
        tuple val(sampleName), path("*_filtered_seurat_sce.rds.gz"), emit: basic_seurat_sce
        tuple val(sampleName), path("*_filtered_sce.rds.gz"), emit: basic_sce, optional: true
        tuple val(sampleName), path("*_unfiltered_seurat_sce.rds.gz"), emit: basic_raw_seurat_sce, optional: true
        tuple val(sampleName), path("*_unfiltered_sce.rds.gz"), emit: basic_raw_sce, optional: true
        tuple val(sampleName), path("*_filtered.h5ad"), emit: basic_anndata, optional: true
        tuple val(sampleName), path("*_unfiltered.h5ad"), emit: basic_raw_anndata, optional: true
        tuple val(sampleName), path("*.pdf"), emit: basic_sce_qc, optional: true

    script:
    outname = flagOn(params.filter, false) ? 'nf-scRNASeq.prefiltered' : 'nf-scRNASeq'
    outprefix = "${sampleName}_${outname}"
    """
    SAMPLES='${sampleName}'
    FEATURES='${featureMatrix}'

    Rscript --vanilla ${sq(rscripts + '/preprocessData.R')} \
       --sample \$SAMPLES \
       --data10X \$FEATURES \
       --baseCond ${params.refName} \
       --annotation ${gtf} \
       --max_mt ${params.maxmtpercent} \
       --min_features ${params.mindetectedfeatures} \
       --hvg_cutoff ${params.hvgcutoff} \
       --output_seurat ${params.output_seurat} \
       --output_sce ${params.output_sce} \
       --output_anndata ${params.output_anndata} \
       --run_sct ${params.run_sct} \
       --out ${outprefix} \
       --marker ${marker} \
       --libpath ${sq(rscripts.toString() + '/')}
    """
}


/************************************************************************
                    STEP 3: Generate overview plots
************************************************************************/
process createOverviewPlots{

    cache 'lenient'
    tag "${sampleName}"

    publishDir "${params.outputDir}/Plots/Overview", mode: 'copy',
    saveAs: {filename ->
        def name = file(filename).getName()
        def lib = name.replaceFirst(/_(nf-scRNASeq)(\.prefiltered)?_.*$/, '')
        "${lib}/${name}"
    }

    input:
        tuple val(sampleName), path(sce), path(rscripts)

    output:
        tuple val(sampleName), path("*overview.pdf"), emit: pdf
        tuple val(sampleName), path("*_QC_*.pdf"), emit: qc_pdf, optional: true
        tuple val(sampleName), path("*.tsv"), emit: tables, optional: true

    script:
    outname = flagOn(params.filter, false) ? "${sampleName}_nf-scRNASeq.prefiltered" : "${sampleName}_nf-scRNASeq"
    """
    Rscript --vanilla ${sq(rscripts + '/create_overview_plots.R')} \
        --sce ${sce} \
        --baseCond ${params.refName} \
        --out ${outname}_overview \
        --format pdf \
        --width 25 \
        --height 10 \
        --libpath ${sq(rscripts.toString() + '/')}
    """
}

/************************************************************************
                    STEP 4: Identify DE Genes
************************************************************************/

process identifyDEGenes{

    cache 'lenient'
    tag "${sampleTag}"

    publishDir "${params.outputDir}", mode: 'copy',
    saveAs: {filename ->
        if (filename.endsWith('.pdf'))       "DE/Plots/${file(filename).getName()}"
        else                                 "DE/GENES/${file(filename).getName()}"
    }

    input:
        tuple val(sampleTag), path(sces), path(rscripts)
        val(check)

    output:
        tuple val(sampleTag), path("*.pdf"), emit: pdf, optional: true
        tuple val(sampleTag), path("*.tsv.gz"), emit: tables, optional: true
        tuple val(sampleTag), path("*.tsv"), emit: summaries, optional: true
        tuple val(sampleTag), path("*.rds.gz"), emit: rds

    script:
    outname = flagOn(params.filter, false) ? "${sampleTag}_nf-scRNASeq.prefiltered" : "${sampleTag}_nf-scRNASeq"
    scelist = (sces instanceof List) ? sces.join(',') : "${sces}"
    """
    Rscript --vanilla ${sq(rscripts + '/identify_de_genes.R')} \
        --sce ${scelist} \
        --baseCond ${params.refName} \
        --format pdf \
        --width 400 \
        --height 300 \
        --organism ${params.organism}\
        --run_pathways ${params.run_pathways} \
        --pcut ${params.pvalcutoff}\
        --fcut ${params.lfccutoff} \
        --out ${outname} \
        --libpath ${sq(rscripts.toString() + '/')}
    """
}


/***********************************************************************
                        MAIN WORKFLOW
************************************************************************/

workflow{
    main:

        // Show the help message if --help is specified
        if (params.help) {
            helpMessage()
            exit 0
        }

        def deOnly = params.de_rds && !(params.de_rds instanceof Boolean)

        def obsoletePipseqOptions = ['pipseq_config_dir', 'pipseq_config_yaml'].findAll { params.containsKey(it) }
        if (obsoletePipseqOptions) {
            error "Options --${obsoletePipseqOptions.join(', --')} were replaced by --pipseq_config"
        }

        if (!deOnly && !params.libraries) {
            log.info("ERROR: --libraries is not specified")
            helpMessage()
            exit 1
        }

        // Check all required input files before they are fed to a channel because this
        // causes the pipeline to fail immediately on AWS before starting up the machines
        if(!deOnly && !file(params.libraries).exists()) {
            log.info("ERROR: the file '${params.libraries}' does not exist")
            exit 2
        }

        // Print execution parameters to stdout
        log.info """
 ======================================================================
 | nf-scRNASeq
 |
 | Version: ${workflow.manifest.version}
 ----------------------------------------------------------------------
 |
 | Mandatory arguments
 |   libraries               : ${params.libraries}
 |
 | Optional arguments
 |   outputDir               : ${params.outputDir}
 |   reportsDir              : ${params.reportsDir}
 |   scriptDirR              : ${params.scriptDirR}
 |   scriptDirPy             : ${params.scriptDirPy}
 |   binDir                  : ${params.binDir}
 |   mapper                  : ${params.mapperbin}
 |   index                   : ${params.mapindex}
 |   reference               : ${params.mapref}
 |   annotation              : ${params.mapanno}
 |   whitelist               : ${params.whitelist}
 |   withQC                  : ${params.runqc}
 |   fastqc_params           : ${params.qcparams}
 |   cellranger_params       : ${params.crparams}
 |   star_params             : ${params.starparams}
 |   idx_params              : ${params.idxparams}
 |   filter                  : ${params.filter}
 |   minReads                : ${params.minReads}
 |   max_mt_percent          : ${params.maxmtpercent}
 |   min_detected_features   : ${params.mindetectedfeatures}
 |   hvg_cutoff              : ${params.hvgcutoff}
 |   output_seurat           : ${params.output_seurat}
 |   output_sce              : ${params.output_sce}
 |   output_anndata          : ${params.output_anndata}
 |   run_sct                 : ${params.run_sct}
 |   pval_cutoff             : ${params.pvalcutoff}
 |   lfc_cutoff              : ${params.lfccutoff}
 |   run_pathways            : ${params.run_pathways}
 |   de_rds                  : ${params.de_rds}
 |   markerfile              : ${params.markerfile}
 |   organism                : ${params.organism}
 |   baseline                : ${params.refName}
 |   pipseq_config           : ${params.pipseqConfig}
 |   pipseq_cb_length        : ${params.pipseqCbLength}
 |   pipseq_umi_length       : ${params.pipseqUmiLength}
 |   pipseq_offset           : ${params.pipseqOffset}
 |   pipseq_umi_offset       : ${params.pipseqUmiOffset}
 |   pipseq_exact            : ${params.pipseqExact}
 |   pipseq_clip3p           : ${params.pipseqClip3p}
 |   pipseq_container        : ${params.pipseqContainer}
 |   r_container             : ${params.rContainer}
 |   stopOnWarnings          : ${params.stopOnWarnings}
 |
 ======================================================================
""".stripIndent()

        if (params.runqc){
            log.info("Running QC")
        }

        if (flagOn(params.filter, false)){
            log.info("Using filtered count output")
        }else{
            log.info("Using raw count output")
        }

        Ch_rscripts = Channel.value(file(params.scriptDirR, checkIfExists: true))

        if (deOnly) {
            def dePaths = params.de_rds.toString().tokenize(',').collect { it.trim() }.findAll { it }
            dePaths.each { path ->
                if (!(path.endsWith('.rds') || path.endsWith('.rds.gz'))) {
                    error "--de_rds must point to .rds or .rds.gz files, got '${path}'"
                }
                if (!file(path).exists()) {
                    error "--de_rds file '${path}' does not exist"
                }
            }
            identifyDEGenes(Channel.value(tuple('de_rds', dePaths.collect { file(it) })).combine(Ch_rscripts), Channel.of([]))
        } else {

        /**********************************************************
                Validate SampleSheet
        ***********************************************************/

        if (file(params.libraries).exists()){
            if (file(params.libraries).isFile()){
                if (params.libraries.endsWith(".csv")){
                    if (file(params.libraries).size() == 0){
                        log.error("SampleSheet is empty!")
                    }
                }else{
                    log.error("SampleSheet is without .csv ending!")
                }
            }else{
                log.error("SampleSheet is not a file!")
            }
        }

        check_samplesheet(Channel.fromPath(params.libraries), Channel.value(file(params.scriptDirPy + 'checkSampleSheet.py', checkIfExists: true)))

        def libraryTypes = file(params.libraries).splitCsv(header: true).collect { canonicalLibraryType(it.LibraryType) }
        def pipseqLibraries = file(params.libraries).splitCsv(header: true).findAll { row -> normalizeChemistry(row.Chemistry) == 'PIPseq' }
        if (pipseqLibraries && params.mapperbin != 'STAR') {
            error("PIPseq CellRanger integration not implemented/verified; use --mapper STAR")
        }
        if (pipseqLibraries.any { canonicalLibraryType(it.LibraryType) == 'CellRanger' }) {
            error("PIPseq libraries must be raw FASTQ libraries with LibraryType GEX")
        }
        if (pipseqLibraries) {
            validate_pipseq_geometry()
        }
        def pipseqConfig = pipseqLibraries ? resolvePipspeakConfig(params.pipseqConfig) : null

        /**********************************************************
                STEP 0: Prepare Input and Indices and run QC
        ***********************************************************/
        
        Ch_csv_raw = check_samplesheet.out.csv.splitCsv(sep: ",", header: true)
            .map { row ->
                row.LibraryType = canonicalLibraryType(row.LibraryType)
                row.Chemistry = normalizeChemistry(row.Chemistry)
                return row
            }

        Ch_pipseq_input = Ch_csv_raw.filter { it.Chemistry == 'PIPseq' }
            .map { row ->
                def key = [row.SampleName, row.Condition, row.Replicate, row.LibraryType].join('|')
                tuple(key, [lane: row.Lane.toString().trim().toInteger(), r1: row.R1, r2: row.R2, sample: row.SampleName, condition: row.Condition, replicate: row.Replicate, librarytype: row.LibraryType, cellnumber: row.CellNumber])
            }
            .groupTuple(by: 0)
            .map { key, entries ->
                def lanes = entries.sort { a, b -> a.lane <=> b.lane }
                def first = lanes[0]
                def id = "${first.sample.replaceAll('_','-')}_${first.condition.replaceAll('_','-')}_${first.replicate.replaceAll('_','-')}"
                def meta = [
                    id: id,
                    sample: first.sample,
                    condition: first.condition,
                    replicate: first.replicate,
                    librarytype: first.librarytype,
                    chemistry: 'PIPseq',
                    cellnumber: first.cellnumber,
                    lanes: lanes.collect { it.lane }
                ]
                tuple(meta, lanes.collect { file(it.r1) }, lanes.collect { file(it.r2) })
            }

        Ch_pipseq_config = pipseqConfig ? Channel.value(tuple(pipseqConfig[0], pipseqConfig[1], pipseqConfig[2])) : Channel.empty()

        preprocessPipseq(Ch_pipseq_input
            .combine(Ch_pipseq_config)
            .combine(Channel.value(file(params.scriptDirPy + 'preparePipseq.py', checkIfExists: true))))

        Ch_pipseq_rows = preprocessPipseq.out.converted.map { meta, convertedR1, convertedR2 ->
            [SampleName: meta.sample, Condition: meta.condition, Replicate: meta.replicate.toString(), Lane: meta.lanes[0].toString(), LibraryType: meta.librarytype, R1: convertedR1.toString(), R2: convertedR2.toString(), CellNumber: meta.cellnumber, Chemistry: 'PIPseq']
        }

        Ch_csv = Ch_csv_raw.filter { it.Chemistry != 'PIPseq' }.mix(Ch_pipseq_rows)

        Ch_csv_GEX_split = Ch_csv.branch{
                raw: (new File(it.R1)).isFile()
                precomputed: (new File(it.R1)).isDirectory()
            }

        Ch_csv_GEX_split.precomputed
            .filter { it.LibraryType != 'CellRanger' }
            .subscribe { log.warn("Sample ${it.SampleName} points R1 at a directory, treating it as precomputed CellRanger output") }

        //Create a list of unique conditions
        Ch_Conditions = Ch_csv.map { row -> row.Condition }.distinct().branch{                 
                multi:  it != params.refName
                single: it == params.refName
        }

        if (params.mapindex){
            if (!file(params.mapindex).exists()){
                if (params.mapperbin == 'CellRanger'){
                    Cellranger_idx(Channel.fromPath(params.mapref), Channel.fromPath(params.mapanno))
                    Ch_mapping_idx = Cellranger_idx.out.idx
                }else if (params.mapperbin == 'STAR'){
                    star_idx(Channel.fromPath(params.mapref), Channel.fromPath(params.mapanno))
                    Ch_mapping_idx = star_idx.out.idxlink
                }
            }else{
                if (params.mapperbin == 'CellRanger'){
                    Ch_mapping_idx = Channel.fromPath(params.mapindex)
                }else if (params.mapperbin == 'STAR'){
                    Ch_mapping_idx = Channel.fromPath(params.mapindex+".idx")
                }
            }
        } else {
            Ch_mapping_idx = Channel.empty()
        }

        if (params.whitelist){
            if (file(params.whitelist).exists()){
                Ch_whitelist = Channel.fromPath(params.whitelist)
            
            } else {         
            create_dummy_whitelist()   
            Ch_whitelist = create_dummy_whitelist.out.dummy
            }
        } else {         
            create_dummy_whitelist()   
            Ch_whitelist = create_dummy_whitelist.out.dummy
        }

        Ch_map_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.CellNumber, row.Chemistry ) }.groupTuple(by: 0)
        
        if (params.mapperbin == 'CellRanger'){
            Ch_map_precomputed = Ch_csv_GEX_split.precomputed.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1)) }.groupTuple(by: 0)            
        }

        if(params.runqc){
            Ch_QC_input = Ch_csv.filter { new File(it.R1).isFile() }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2)) }.groupTuple(by: 0)
            qc_raw(Ch_QC_input)
            mqc(qc_raw.out.zip.collect())
        }

        /**********************************************************
                STEP 0.1: Concat Lanes
        ***********************************************************/

        concat_lanes(Ch_map_input)

        /**********************************************************
                STEP 1: Count Reads
        ***********************************************************/

        if (params.mapperbin == 'CellRanger'){
            runCellrangerCount(concat_lanes.out.combine( Ch_mapping_idx ))
            useCellrangerData(Ch_map_precomputed)
            if (flagOn(params.filter, false)){
                Ch_cell_data = runCellrangerCount.out.cell_data_filtered.mix(useCellrangerData.out.cell_data_from_precomputed_filtered)
            } else {
                Ch_cell_data = runCellrangerCount.out.cell_data_raw.mix(useCellrangerData.out.cell_data_from_precomputed_raw)
            }
        }else if (params.mapperbin == 'STAR'){
            star_mapping(concat_lanes.out.combine( Ch_mapping_idx ).combine( Ch_whitelist ))                    
            if (flagOn(params.filter, false)){
                Ch_cell_data = star_mapping.out.cell_data_filtered                
            } else {
                Ch_cell_data = star_mapping.out.cell_data_raw
            }
        }
    
        /**************************************************************
            STEP 2: Generate SingleCellExperiment object
        ***************************************************************/

        Ch_gtf = Channel.fromPath(params.mapanno)

        Ch_preprocess_input = Ch_cell_data.combine(Ch_gtf)
                        .combine(Channel.value(file(params.markerfile, checkIfExists: true)))
                        .combine(Ch_rscripts)
                        .map { sampleName, featureMatrix, gtf, marker, rscripts ->
                        tuple(sampleName, featureMatrix, gtf, marker, rscripts)
                        }

        preprocessSingleCellData(Ch_preprocess_input)

        /**************************************************************
            STEP 3: Generate overview plots
        ***************************************************************/
        createOverviewPlots(preprocessSingleCellData.out.basic_seurat_sce.combine(Ch_rscripts))
        
        /**************************************************************
            STEP 4: Run DE Analysis for Genes
        ***************************************************************/
        Ch_multi_conditions = Ch_Conditions.multi.collect().ifEmpty([])

        Ch_de_input = preprocessSingleCellData.out.basic_seurat_sce
            .map { sampleName, sce -> tuple(sampleName, sce) }
            .toList()
            .filter { rows -> rows.size() > 0 }
            .map { rows ->
                def sorted = rows.sort { row -> row[1].getName() }
                tuple('allSamples', sorted.collect { row -> row[1] })
            }
            .combine(Ch_rscripts)

        identifyDEGenes(Ch_de_input, Ch_multi_conditions)
        }
}
