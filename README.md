# RNAseq-variant-calling-nf
`RNAseq-variant_calling-nf` is a A Nextflow workflow to call variants with GATK from bulk RNAseq data  

# Usage
`nextflow run RNAseq-variant-calling-nf --fastq [path to fastqs] --project [project ID]`

# Required arguments:
+ `--fastq` 
    + Path to fastqs, default is `./Fastq`
+ `--project`
    + Project name. E.g., `'19482R'`. Will be given to final vcf file

# Outline
This workflow will:
+ Make channel from paired end reads
+ Remove optical duplicates with Clumpify
+ Trim adapters with Cutadapt (hardcoded adapter sequence)
+ Map with STAR. BAMs published to outdir/bams
+ Index BAMS and calculate idxstats with Samtools. Also published to oudir/bams
+ featureCounts. Published to outdir/feature_counts
+ RSEM gene and isoform counts are published to outdir/rsem
+ RNAseq metrics with Picard. Published to outdir/rnaseq_metrics
+ Fliter alignments for unmapped, secondary, qcfail, and supplementary        
+ Mark duplicates with Picard and reindex.
+ Generate read coverage with USeq and convert to BigWig. Published to outdir/coverage 
+ Split and trim intron junctions
+ Indel realigner
+ Base recalibration
+ Pad BED file
+ Haplotype variant calls. Published to outdir/haplotype_vcf
+ Variant filtering
+ Merge raw variants
+ Generate genomic VCF for joint genotyping
+ Joint genotyping. Published to outdir/final

# Results
Workflow results are published to `--outdir`, whose default is `./results`  
These include:
+ BAMs from STAR to `./results/bams`
+ Samtools idxstats to `./results/bams`
+ featureCounts to `./results/feature_counts`
+ RSEM to `./results/rsem`
+ RNAseq metrics to `./results/rnaseq_metrics`
+ BigWigs to `./results/coverage`
+ Per sample GATK haplotype calls to `./results/haplotype_vcf`
+ Joint genotyping with GATK to `./results/final`

# Defaults
Default values for paths to applications and reference files are in `nextflow.config`.  
Default organism is human with Ensembl v96 annotation
