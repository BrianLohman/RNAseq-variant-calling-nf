#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info"""
    -----------------------------------------------------------------------------
    RNAseq_variants: a workflow for calling variants from RNAseq data
    =============================================================================

    Required arguments:
    -------------------

    --reads         Full path to directory with reads. Default: ./fastq

    --outdir        Path to publish results. Default: ./results

    --genome        Full path to STAR reference. Default: Hg38 in core group space

    --rsem_index    Full path to RSEM index. Default: Hg38 in core group space

    
    Description:
    ------------
    Index BAMS (BAMS and indicies are published to params.outdir/bams)

    Run RSEM to get gene and isoform counts (saved to params.outdir/rsem_out)

    -----------------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// Set default required params
params.reads = "./fastq/*_R{1,2}_001.fastq.gz"
params.outdir = './results'
params.genome='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/star100'
params.rsem_index='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/rsem/RSEM'
params.gtf='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/Homo_sapiens.GRCh38.96.gtf'
params.refflat='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/Homo_sapiens.GRCh38.96.refflat'
params.riboint='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/Homo_sapiens.GRCh38.96.rRNA.interval'
params.chromsize='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/chrom.sizes'
params.ref='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.goldindels='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.hiconfsnps='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/GATK/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
params.dbsnp='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/GATK/dbsnp_144.hg38.vcf.gz'
params.exons='/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/Homo_sapiens.GRCh38.96.mergedExons.bed'


// Logging
log.info("\n")
log.info("Reads directory (--reads)         :${params.reads}")
log.info("Genome          (--genome)        :${params.genome}")
log.info("RSEM index      (--rsem_index)    :${params.rsem_index}")
log.info("Results         (--outdir)        :${params.outdir}")
log.info("GTF             (--gft)           :${params.gtf}")
log.info("REFFLAT         (--refflat)       :${params.refflat}")
log.info("RIBOINT         (--riboint)       :${params.riboint}")
log.info("CHROMSIZE       (--chromsize)     :${params.chromsize}")
log.info("REF (fasta)     (--ref)           :${params.ref}")
log.info("GOLDINDESL      (--goldindels)    :${params.goldindels}")
log.info("HIGHCONFSNPS    (--highconfsnps)  :${params.hiconfsnps}")
log.info("DBSNP           (--dbsnp)         :${params.dbsnp}")
log.info("EXONS           (--exons)         :${params.exons}")

// Read pair channel: tuple of pair_id and paired fastqs
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

// Remove optical duplicates from fastq
process dedup {
    tag "$pair_id"

    input:
      tuple val(pair_id), path(reads)

    output:
      tuple val(pair_id), path("${pair_id}.read*.fastq.gz"), emit: deduped_reads
      path("${pair_id}.clumpify.out.txt")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/BBmap/v38.34/clumpify.sh \
      in1=${pair_id}_R1_001.fastq.gz \
      in2=${pair_id}_R2_001.fastq.gz \
      out1=${pair_id}.read1.fastq.gz \
      out2=${pair_id}.read2.fastq.gz \
      dupedist=10000 \
      dedupe=t \
      optical=t 2> ${pair_id}.clumpify.out.txt
      """
}

// Trim adapters
process trim {
    module 'cutadapt/2.10'
    tag "$pair_id"

    input:
      tuple val(pair_id), path(deduped_reads)

    output:
      tuple val(pair_id), path("${pair_id}.trim*.fastq.gz"), emit: trimmed_reads
      path("${pair_id}.cutadapt.out.txt")

    script:
      """
      cutadapt \
      -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
      -o ${pair_id}.trim1.fastq.gz \
      -p ${pair_id}.trim2.fastq.gz \
      -j 8 \
      -m 25 \
      ${pair_id}.read1.fastq.gz ${pair_id}.read2.fastq.gz > ${pair_id}.cutadapt.out.txt
      """
}

// Maps each read-pair by using STAR
process star {
    tag "$pair_id"
      
    input:
      path(genome)
      tuple val(pair_id), path(reads)
  
    output:
      tuple val(pair_id), path("${pair_id}.bam"), emit: bam
      tuple val(pair_id), path("${pair_id}_aligned_to_transcriptome.out.bam"), emit: rsem_input

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/STAR/2.7.0f/STAR --genomeDir $genome \
      --readFilesIn $reads \
      --readFilesCommand zcat \
      --runThreadN 12 \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 34359720776 \
      --outBAMsortingBinsN 100 \
      --quantMode TranscriptomeSAM \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outSAMstrandField intronMotif \
      --outSAMattributes NH HI NM MD
      mv Aligned.sortedByCoord.out.bam ${pair_id}.bam
      mv Aligned.toTranscriptome.out.bam ${pair_id}_aligned_to_transcriptome.out.bam
      """
}

// Index BAMs
process index {
    tag "${pair_id}"
    publishDir "${params.outdir}/bams", mode:"copy"

    input:
      tuple val(pair_id), path(bam)

    output:
      tuple val(pair_id), path("${pair_id}.bam"), path("${pair_id}.bam.bai"), emit: indexed_bam
      path("${pair_id}.idxstats")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index ${pair_id}.bam
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools idxstats ${pair_id}.bam | sort -V > ${pair_id}.idxstats
      """
}

// FeatureCounts
process featureCounts {
    tag "${pair_id}"

    publishDir "${params.outdir}/feature_counts", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)

    output:
      path("${pair_id}.counts.gz")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/Subread/1.5.1/bin/featureCounts -T 8 -p -C -s 2 --largestOverlap -a /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release96/Homo_sapiens.GRCh38.96.gtf -o ${pair_id}.counts ${pair_id}.bam
      gzip ${pair_id}.counts
      """
}

// RSEM to get gene and isoform counts
process rsem {
    tag "$pair_id"
    publishDir "${params.outdir}/rsem_out", mode:"copy"
 
    input:
      path(params.rsem_index)
      tuple val(pair_id), path(rsem_bam)
 
    output:
      path("${pair_id}.genes.results")
      path("${pair_id}.isoforms.results")
 
    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/rsem/1.3.1/rsem-calculate-expression \
      --alignments --paired-end --num-threads 8 \
      --no-bam-output --quiet \
      ${rsem_bam} ${params.rsem_index} ${pair_id}
      """
}

workflow {
    dedup(read_pairs_ch)
    trim(dedup.out.deduped_reads)
    star(params.genome, trim.out.trimmed_reads)
    index(star.out.bam)
    featureCounts(index.out.indexed_bam)
    rsem(params.rsem_index, star.out.rsem_input)
}
