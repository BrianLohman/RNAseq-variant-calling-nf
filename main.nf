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
    Call variants with GATK from bulk RNAseq data. Defaults to human Ensembl v96.
    All reference files can be set via params in nextflow.config.

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
        Generate read coverage with USeq and convert to BigWig. Published to outdir/coverage 
        Split and trim intron junctions
        Indel realigner
        Base recalibration
        Pad BED file
        Haplotype variant calls. Published to outdir/haplotype_vcf
        Variant filtering
        Merge raw variants
        Generate genomic VCF for joint genotyping
        Joint genotyping. Published to outdir/final

    -----------------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// Check for required user input
if ( !params.project ) { exit 1, "--project is not defined" }

// Logging
log.info("\n")
log.info("Fastq directory       (--fastq)                 :${params.fastq}")
log.info("Genome                (--genome)                :${params.genome}")
log.info("RSEM index            (--rsem_index)            :${params.rsem_index}")
log.info("Results               (--outdir)                :${params.outdir}")
log.info("GTF                   (--gft)                   :${params.gtf}")
log.info("REFFLAT               (--refflat)               :${params.refflat}")
log.info("RIBOINT               (--riboint)               :${params.riboint}")
log.info("CHROMSIZE             (--chromsize)             :${params.chromsize}")
log.info("REF (fasta)           (--ref)                   :${params.ref}")
log.info("REF INDEX             (--ref_index)             :${params.ref_index}")
log.info("GOLDINDESL            (--goldindels)            :${params.goldindels}")
log.info("HIGHCONFSNPS          (--highconfsnps)          :${params.hiconfsnps}")
log.info("HIGHCONFSNPS INDEX    (--highconfsnps_index)    :${params.hiconfsnps_index}")
log.info("DBSNP                 (--dbsnp)                 :${params.dbsnp}")
log.info("DBSNP INDEX           (--dbsnp_index)           :${params.dbsnp_index}")
log.info("EXONS                 (--exons)                 :${params.exons}")
log.info("DICT                  (--dict)                  :${params.dict}")
log.info("HAPMAP                (--hapmap)                :${params.hapmap}")
log.info("HAPMAP INDEX          (--hapmap_index)          :${params.hapmap_index}")
log.info("OMNI                  (--omni)                  :${params.omni}")
log.info("OMNI INDEX            (--omni_index)            :${params.omni_index}")
log.info("PROJECT               (--project)               :${params.project}")


// Read pair channel: tuple of pair_id and paired fastqs
Channel
    .fromFilePairs( params.fastq )
    .ifEmpty { error "Cannot find any reads matching: ${params.fastq}" }
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

// Maps each read-pair with STAR
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
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/STAR/2.7.0f/STAR --runMode alignReads \
      --genomeDir $genome \
      --twopassMode Basic \
      --readFilesIn $reads \
      --readFilesCommand zcat \
      --runThreadN ${task.cpus} \
      --alignIntronMax 100000 \
      --alignIntronMin 20 \
      --alignMatesGapMax 100000 \
      --outFilterMismatchNmax 20 \
      --outFilterMismatchNoverLmax 0.3 \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 25000000000 \
      --outSAMmapqUnique 60 \
      --outSAMattrRGline ID:${pair_id} LB:${pair_id} PL:Illumina PU:${pair_id} SM:${pair_id} CN:HCI \
      --quantMode TranscriptomeSAM \
      --quantTranscriptomeBan IndelSoftclipSingleend \
      --outWigType None 
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

// RSEM: gene and isoform counts
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
// Extract regions from BED file if specified in params.bed
process filter_alignments {
    tag "${pair_id}"

    input:
      tuple val(pair_id), path(bam), path(bai)
      path(bed)

    output:
      tuple val(pair_id), path("${pair_id}.bam"), path("${pair_id}.bam.bai"), emit: filtered_bam

    script:
      if( params.bed == 'NO_FILE' ) {
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools view -L $bed -o ${pair_id}_subset.bam $bam
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools view -F 2820 -@ 8 -b -o ${pair_id}.filter.bam ${pair_id}_subset.bam
      """
      } else {
      """
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index -b -@ 8 ${pair_id}.filter.bam ${pair_id}.filter.bai
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools view -F 2820 -@ 8 -b -o ${pair_id}.filter.bam $bam
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index -b -@ 8 ${pair_id}.filter.bam ${pair_id}.filter.bai
      """
      }
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

// Generate read coverage and passing BED file
process useq {
    tag "${pair_id}"

    stageInMode 'copy'
    
    publishDir "${params.outdir}/coverage", mode:"copy"

    input:
      tuple val(pair_id), path(bam), path(bai)
      path(exons)

    output:
      path("${pair_id}.region_stats.txt.gz")
      path("${pair_id}.stats.json.gz")
      path("${pair_id}.stats.txt")
      path("*.bw")
      tuple val(pair_id), path("${pair_id}_Pass.bed.gz"), emit: passing_bed

    script:
      """
      java -Xmx${task.memory.giga}g -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/useq/9.1.3/Apps/Sam2USeq \
      -f $bam \
      -v H_sapiens_Dec_2013 \
      -m 10 -a 1000 -r \
      -b $exons \
      -x 1000 -c 20 \
      -p ${pair_id}.region_stats.txt.gz \
      -j ${pair_id}.stats.json.gz \
      -o ${pair_id}.stats.txt \
      -n "${pair_id}"

      java -Xmx${task.memory.giga}g -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/useq/9.1.3/Apps/USeq2UCSCBig \
      -u ./ \
      -d /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/UCSC
      """
}

// Split and trim intron junctions
process intron_junctions {
    tag "${pair_id}"

    input:
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(pair_id), path(bam), path(bai)

    output:
      tuple val(pair_id), path("${pair_id}.split.bam"), path("${pair_id}.split.bam.bai"), emit: split_bam

    script:
      """
      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T SplitNCigarReads \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -I $bam \
      -o ${pair_id}.split.bam \
      -U ALLOW_N_CIGAR_READS \
      -drf DuplicateRead 

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index ${pair_id}.split.bam
      """
}

// Indel realigner
process indel_realigner {
    tag "${pair_id}"

    input:
      path(ref)
      path(ref_index)
      path(dict)
      path(goldindels)
      path(goldindels_index)
      tuple val(pair_id), path(bam), path(bai)

    output:
      tuple val(pair_id), path("${pair_id}.realign.bam"), path("${pair_id}.realign.bam.bai"), emit: realigned_bam

    script:
      """
      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -I $bam \
      -o ${pair_id}.intervals \
      --known $goldindels \
      -nt ${task.cpus}

      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T IndelRealigner \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -targetIntervals ${pair_id}.intervals \
      -I $bam \
      -o ${pair_id}.realign.bam \
      -drf DuplicateRead \
      -known $goldindels 

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index ${pair_id}.realign.bam
      """
}

// Base recalibration
process base_recalibration {
    tag "${pair_id}"

    input:
      path(ref)
      path(ref_index)
      path(dict)
      path(hiconfsnps)
      path(hiconfsnps_index) 
      path(goldindels)
      path(goldindels_index)
      tuple val(pair_id), path(bam), path(bai)

    output:
      tuple val(pair_id), path("${pair_id}.final.bam"), path("${pair_id}.final.bam.bai"), emit: final_bam

    script:
      """
      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T BaseRecalibrator \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -knownSites $hiconfsnps \
      -knownSites $goldindels \
      -I $bam \
      -o ${pair_id}.realign.grp \
      -nct ${task.cpus}

      java -Xmx${task.memory.giga}g -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T PrintReads \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -I $bam \
      -BQSR ${pair_id}.realign.grp \
      -o ${pair_id}.final.bam \
      -drf DuplicateRead \
      -nct ${task.cpus} 

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/samtools/1.8/samtools index ${pair_id}.final.bam
      """
}

// Haplotype variant calling
process variant_calling {
    tag "${pair_id}"

    publishDir "${params.outdir}/haplotype_vcf", mode:"copy"

    input:
      path(ref)
      path(ref_index)
      path(dbsnp)
      path(dbsnp_index) 
      path(dict)
      path(chromsize)
      tuple val(pair_id), path(passing_bed)
      tuple val(pair_id), path(bam), path(bai)
      
    output:
      tuple val(pair_id), path("${pair_id}.vcf.gz"), path("${pair_id}.vcf.gz.tbi"), emit: vcf

    script:
      """
      # split bed file for manual parallel processing
      zcat $passing_bed | \
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/bedtools/2.22.1/bedtools slop -g $chromsize -b 20 -i - | \
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/bedtools/2.22.1/bedtools sort -i - > "${pair_id}_Pass.ext20.bed"
      split -a 1 -n l/12 --additional-suffix=.bed "${pair_id}_Pass.ext20.bed" "${pair_id}_call."
      gzip "${pair_id}_Pass.ext20.bed"
      
      # generate variant calls
      parallel -k --delay 10 \
      java -Xmx4G -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T HaplotypeCaller \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      --dbsnp $dbsnp \
      -I $bam \
      -dontUseSoftClippedBases \
      -stand_call_conf 20.0 \
      -mmq 1 \
      --min_base_quality_score 20 \
      -L {} \
      -o {.}.raw.vcf.gz \
      ':::' ${pair_id}_call.*.bed

      # variant filtering
      parallel -k --delay 10 \
      java -Xmx4G -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T VariantFiltration \
      -R $ref \
      -V {.}.raw.vcf.gz \
      -window 35 \
      -cluster 3 \
      -filterName FS \
      -filter "'FS > 30.0'" \
      -filterName QD \
      -filter "'QD < 2.0'" \
      -o {.}.vcf.gz \
      ':::' ${pair_id}_call.*.bed

      # merge raw variants
      java -Xmx20G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar SortVcf \
      SD=$dict \
      O=${pair_id}.vcf.gz \
      I=${pair_id}_call.a.vcf.gz \
      I=${pair_id}_call.b.vcf.gz \
      I=${pair_id}_call.c.vcf.gz \
      I=${pair_id}_call.d.vcf.gz \
      I=${pair_id}_call.e.vcf.gz \
      I=${pair_id}_call.f.vcf.gz \
      I=${pair_id}_call.g.vcf.gz \
      I=${pair_id}_call.h.vcf.gz \
      I=${pair_id}_call.i.vcf.gz \
      I=${pair_id}_call.j.vcf.gz \
      I=${pair_id}_call.k.vcf.gz \
      I=${pair_id}_call.l.vcf.gz \

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${pair_id}.vcf.gz
      """
}

// Generate genomic VCF for joint genotyping
process genomic_vcf {
    tag "${pair_id}"

    publishDir "${params.outdir}/genomic_vcf", mode:"copy"
    
    input:
      path(exons)
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(pair_id), path(bam), path(bai)

    output:
      path("${pair_id}.g.vcf.gz"), emit: gvcf
      path("${pair_id}.g.vcf.gz.tbi"), emit: gvcf_index

    script:
      """
      # split file by lines into chunks
      split -a 1 -n l/12 --additional-suffix=.bed $exons exon_call.

      # haplotype call
      parallel -k --delay 10 \
      java -Xmx4G -XX:ParallelGCThreads=${task.cpus} \
      -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar -T HaplotypeCaller \
      -jdk_deflater -jdk_inflater \
      -R $ref \
      -I $bam \
      -dontUseSoftClippedBases \
      --genotyping_mode DISCOVERY \
      --emitRefConfidence GVCF \
      -stand_call_conf 20.0 \
      -mmq 1 \
      --min_base_quality_score 20 \
      -L {} \
      --interval_padding 20 \
      -o {.}.g.vcf.gz \
      ':::' exon_call.*.bed

      # Concatenate gVCF
      java -Xmx20G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar SortVcf \
      SD=$dict \
      O=${pair_id}.g.vcf.gz \
      I=exon_call.a.g.vcf.gz \
      I=exon_call.b.g.vcf.gz \
      I=exon_call.c.g.vcf.gz \
      I=exon_call.d.g.vcf.gz \
      I=exon_call.e.g.vcf.gz \
      I=exon_call.f.g.vcf.gz \
      I=exon_call.g.g.vcf.gz \
      I=exon_call.h.g.vcf.gz \
      I=exon_call.i.g.vcf.gz \
      I=exon_call.j.g.vcf.gz \
      I=exon_call.k.g.vcf.gz \
      I=exon_call.l.g.vcf.gz \

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${pair_id}.g.vcf.gz
      """
}

// Merge sample gVCF files
process joint_genotype {

    input:
      path(exons)
      path(ref)
      path(ref_index)
      path(dict)
      path(dbsnp)
      path(dbsnp_index)
      val(project)
      path(all_gvcf)
      path(all_gvcf_index)

    output:
      tuple val(project), path("${project}.raw.vcf.gz"), path("${project}.raw.vcf.gz.tbi"), emit: raw_merged_vcf

    script:
      """
      # split exon file by lines into chunks
      split -a 1 -n l/16 --additional-suffix=.bed $exons exon_call.

      # Merge sample gVCF file
      parallel -v --delay 30 \
      java -Xmx6G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T CombineGVCFs \
      -R $ref \
      -L exon_call.{}.bed \
      --interval_padding 20 \
      --variant ${all_gvcf.join(' --variant ')} \
      -o ${project}.{}.g.vcf.gz \
      ':::' a b c d e f g h i j k l m n o p

      # genotype
      parallel -k --delay 30 \
      java -Xmx6G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T GenotypeGVCFs \
      -R $ref \
      --dbsnp $dbsnp \
      --variant ${project}.{}.g.vcf.gz \
      -o ${project}.{}.raw.vcf.gz \
      ':::' a b c d e f g h i j k l m n o p

      # Combine chunks
      java -Xmx20G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/picard/2.9.0/picard.jar SortVcf \
      SD=$dict \
      O=${project}.raw.vcf.gz \
      I=${project}.a.raw.vcf.gz \
      I=${project}.b.raw.vcf.gz \
      I=${project}.c.raw.vcf.gz \
      I=${project}.d.raw.vcf.gz \
      I=${project}.e.raw.vcf.gz \
      I=${project}.f.raw.vcf.gz \
      I=${project}.g.raw.vcf.gz \
      I=${project}.h.raw.vcf.gz \
      I=${project}.i.raw.vcf.gz \
      I=${project}.j.raw.vcf.gz \
      I=${project}.k.raw.vcf.gz \
      I=${project}.l.raw.vcf.gz \
      I=${project}.m.raw.vcf.gz \
      I=${project}.n.raw.vcf.gz \
      I=${project}.o.raw.vcf.gz \
      I=${project}.p.raw.vcf.gz

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${project}.raw.vcf.gz
      """
}

// Hard filter excess heterozygosity
// phred score of 54.69 corresponds to a z-score -4.5
process filter_het {
    
    input:
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(project), path(raw_merged_vcf), path(raw_merged_vcf_index)

    output:
      tuple val(project), path("${project}.excessHet.vcf.gz"), path("${project}.excessHet.vcf.gz.tbi"), emit: excessive_het

    script:
      """
      java -Xmx24G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T VariantFiltration \
      -R $ref \
      -V ${raw_merged_vcf} \
      -filter "ExcessHet > 54.69" \
      -filterName ExcessHet \
      -o ${project}.excessHet.vcf.gz
      
      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${project}.excessHet.vcf.gz
      """
}

// SNP Variant Recalibration
process snp_recalibration {
    stageInMode 'copy'

    input:
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(project), path(excessive_het), path(excessive_het_index)
      path(hapmap)
      path(hapmap_index)
      path(omni)
      path(omni_index)
      path(hiconfsnps)
      path(hiconfsnsp_index)
      path(dbsnp)
      path(dbsnp_index)

    output:
      path("${project}.snp.tranches"), emit: snp_tranches
      path("${project}.snp.recal"), emit: snp_recal
      path("${project}.snp.plots.R")

    script:
      """
      java -Xmx24G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T VariantRecalibrator \
      -R $ref \
      -input $excessive_het \
      -mode SNP \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 $hapmap \
      -resource:omni,known=false,training=true,truth=true,prior=12 $omni \
      -resource:1000G,known=false,training=true,truth=false,prior=10 $hiconfsnps \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 $dbsnp \
      -an QD \
      -an ReadPosRankSum \
      -an FS \
      -an SOR \
      --trustAllPolymorphic \
      --maxGaussians 6 \
      -tranchesFile ${project}.snp.tranches \
      -recalFile ${project}.snp.recal \
      -rscriptFile ${project}.snp.plots.R
      """
}

// INDEL Variant Recalibration
process indel_recalibration {

    input:
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(project), path(excessive_het), path(excessive_het_index)
      path(goldindels)
      path(goldindels_index)
      path(dbsnp)
      path(dbsnp_index)

    output:
      path("${project}.indel.tranches"), emit: indel_tranches
      path("${project}.indel.recal"), emit: indel_recal
      path("${project}.indel.plots.R")

    script:
      """
      java -Xmx24G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T VariantRecalibrator \
      -R $ref \
      -input $excessive_het \
      -mode INDEL \
      -resource:mills,known=false,training=true,truth=true,prior=12 $goldindels \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 $dbsnp \
      -an QD \
      -an ReadPosRankSum \
      -an FS \
      -an SOR \
      --trustAllPolymorphic \
      --maxGaussians 4 \
      -tranchesFile ${project}.indel.tranches \
      -recalFile ${project}.indel.recal \
      -rscriptFile ${project}.indel.plots.R
      """
}

// Apply SNP Recalibration
process apply_snp_recalibration {

    input:
      path(ref)
      path(ref_index)
      path(dict)
      tuple val(project), path(excessive_het), path(excessive_het_index)
      path(snp_tranches)
      path(snp_recal)

    output:
      tuple val(project), path("${project}.snp.recal.vcf.gz"), path("${project}.snp.recal.vcf.gz.tbi"), emit: snp_recal_vcf

    script:
      """
      java -Xmx24G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T ApplyRecalibration \
      -R $ref \
      -input $excessive_het \
      -mode SNP \
      -tranchesFile $snp_tranches \
      -recalFile $snp_recal \
      --ts_filter_level 95 \
      -o ${project}.snp.recal.vcf.gz

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${project}.snp.recal.vcf.gz
      """
}

// Apply INDEL Recalibration
process apply_indel_recalibration {

    publishDir "${params.outdir}/final", mode:"copy"

    input:
      path(ref)
      path(ref_index)
      path(dict)
      val(project)
      tuple val(project), path(snp_recal_vcf), path(snp_recal_vcf_index)
      path(indel_tranches)
      path(indel_recal)

    output:
      path("${project}.vqsr.vcf.gz")
      path("${project}.vqsr.vcf.gz.tbi")

    script:
      """
      java -Xmx24G -jar /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/gatk/3.8/GenomeAnalysisTK.jar \
      -jdk_deflater -jdk_inflater \
      -T ApplyRecalibration \
      -R $ref \
      -input $snp_recal_vcf \
      -mode INDEL \
      -tranchesFile $indel_tranches \
      -recalFile $indel_recal \
      --ts_filter_level 95 \
      -o ${project}.vqsr.vcf.gz

      /uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/app/tabix/1.8/tabix -fp vcf ${project}.vqsr.vcf.gz
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
    filter_alignments(index.out.indexed_bam, params.bed)
    mark_duplicates(filter_alignments.out.filtered_bam)
    useq(mark_duplicates.out.mkdup_bam, params.exons)
    intron_junctions(params.ref, params.ref_index, params.dict, mark_duplicates.out.mkdup_bam)
    indel_realigner(params.ref, params.ref_index, params.dict, params.goldindels, params.goldindels_index, intron_junctions.out.split_bam)
    base_recalibration(params.ref, params.ref_index, params.dict, params.hiconfsnps, params.hiconfsnps_index, params.goldindels, params.goldindels_index, indel_realigner.out.realigned_bam)
    variant_calling(params.ref, params.ref_index, params.dbsnp, params.dbsnp_index, params.dict, params.chromsize, useq.out.passing_bed, base_recalibration.out.final_bam)
    genomic_vcf(params.exons, params.ref, params.ref_index, params.dict, base_recalibration.out.final_bam)
    joint_genotype(params.exons, params.ref, params.ref_index, params.dict, params.dbsnp, params.dbsnp_index, params.project, genomic_vcf.out.gvcf.collect(), genomic_vcf.out.gvcf_index.collect())
    filter_het(params.ref, params.ref_index, params.dict, joint_genotype.out.raw_merged_vcf)
    snp_recalibration(params.ref, params.ref_index, params.dict, filter_het.out.excessive_het, params.hapmap, params.hapmap_index, params.omni, params.omni_index, params.hiconfsnps, params.hiconfsnps_index, params.dbsnp, params.dbsnp_index)
    indel_recalibration(params.ref, params.ref_index, params.dict, filter_het.out.excessive_het, params.goldindels, params.goldindels_index, params.dbsnp, params.dbsnp_index)
    apply_snp_recalibration(params.ref, params.ref_index, params.dict, filter_het.out.excessive_het, snp_recalibration.out.snp_tranches, snp_recalibration.out.snp_recal)
    apply_indel_recalibration(params.ref, params.ref_index, params.dict, params.project, apply_snp_recalibration.out.snp_recal_vcf, indel_recalibration.out.indel_tranches, indel_recalibration.out.indel_recal)
}
