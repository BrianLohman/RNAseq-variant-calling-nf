#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info"""
    -----------------------------------------------------------------------------
    RNAseq-variant-calling: a workflow for calling variants from RNAseq data
    =============================================================================

    Required arguments:
    -------------------

    --fastq         Full path to directory with reads. Default: ./fastq

    --outdir        Path to publish results. Default: ./results

    
    Description:
    ------------
    Call variants with GATK from bulk RNAseq data. Defaults to human Ensembl v96 but all reference files can be set via params.
    Outline:
        Make channel from paired end reads
        Remove optical duplicates with Clumpify
        Trim adapters with Cutadapt (hardcoded adapter sequence)
        Map with STAR. BAMs published to outdir/bams
        Index BAMS and calculate idxstats with Samtools. Also published to oudir/bams
        featureCounts. Published to outdir/feature_counts
        RSEM gene and isoform counts are published to outdir/rsem_out
        RNAseq metrics with Picard. Published to outdir/rnaseq_metrics
        Fliter alignments for unmapped, secondary, qcfail, and supplementary        
        Mark duplicates with Picard and reindex. TODO: currently indexing with samtools because Picard not writing bai
        
        TODO:
            Generate read coverage and passing bed file
            Convert USeq to BigWig
            Split and trim intron junctions
            Indel realigner
            Base recalibration
            Pad BED file
            Haplotype variant calls
            Variant filtering
            Merge raw variants
            Generate genomic VCF for joint genotyping
            Joint genotyping

            Add option to extract chroms or regions after STAR?

    -----------------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// Set default required params
params.fastq = "./fastq/*_R{1,2}_001.fastq.gz"
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
log.info("Fastq directory (--fastq)         :${params.fastq}")
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
    .fromFilePairs( params.fastq )
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
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index $bam
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools idxstats ${pair_id}.bam | sort -V > ${pair_id}.idxstats
      """
}

// FeatureCounts
process feature_counts {
    tag "${pair_id}"
    publishDir "${params.outdir}/feature_counts", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)
      path(gtf)

    output:
      path("${pair_id}.counts.gz")

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/Subread/1.5.1/bin/featureCounts \
      -T 8 \
      -p \
      -C \
      -s 2 \
      --largestOverlap \
      -a $gtf \
      -o ${pair_id}.counts \
      $bam
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
      --alignments \
      --paired-end \
      --num-threads 8 \
      --no-bam-output \
      --quiet \
      ${rsem_bam} \
      ${params.rsem_index} \
      ${pair_id}
      """
}

// RNAseq metrics
process rnaseq_metrics {
    tag "${pair_id}"

    publishDir "${params.outdir}/rnaseq_metrics", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)
      path(refflat)
      path(riboint)

    output:
      path("${pair_id}.rna_metrics")

    script:
      """
      java -Xmx${task.memory.giga}g -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar \
      CollectRnaSeqMetrics \
      REF_FLAT=$refflat \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$riboint \
      I=$bam \
      O=${pair_id}.rna_metrics
      """
}

// Filter alignments and reindex
// Unmapped, secondary, qcfail, supplementary = 2820 = 0xb04
process filter_alignments {
    tag "${pair_id}"

    input:
      tuple val(pair_id), path(bam), path(bai)

    output:
      tuple val(pair_id), path("${pair_id}.bam"), path("${pair_id}.bam.bai"), emit: filtered_bam

    script:
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools view -F 2820 -@ 8 -b -o ${pair_id}.filter.bam $bam
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index -b -@ 8 ${pair_id}.filter.bam ${pair_id}.filter.bai
      """
}   

// Mark duplicates and reindex
process mark_duplicates {
    tag "${pair_id}"

    input:
      tuple val(pair_id), path(bam), path(bai)

    output:
      tuple val(pair_id), path("${pair_id}.mkdup.bam"), path("${pair_id}.mkdup.bam.bai"), emit: mkdup_bam

    script:
      """
      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar MarkDuplicates \
      I=$bam \
      O=${pair_id}.mkdup.bam \
      M=${pair_id}.markduplicates.txt \
      REMOVE_DUPLICATES=false
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index ${pair_id}.mkdup.bam
      """
}

workflow {
    dedup(read_pairs_ch)
    trim(dedup.out.deduped_reads)
    star(params.genome, trim.out.trimmed_reads)
    index(star.out.bam)
    feature_counts(index.out.indexed_bam, params.gtf)
    rsem(params.rsem_index, star.out.rsem_input)
    rnaseq_metrics(index.out.indexed_bam, params.refflat, params.riboint)
    filter_alignments(index.out.indexed_bam)
    mark_duplicates(filter_alignments.out.filtered_bam)
}
